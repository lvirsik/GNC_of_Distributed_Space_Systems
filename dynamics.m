classdef dynamics
    methods (Static)
        function statedot = two_body_dynamics(t, state, simulation_settings)
            % Take in state vector and return statedot vector (accelerations)
            r = state(1:3);
            v = state(4:6);

            a = -constants.mu * r / (norm(r) ^ 3);

            if simulation_settings.J2
                a_j2_x = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (1 - (5 * (r(3) / norm(r)) ^ 2)) * r(1);
                a_j2_y = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (1 - (5 * (r(3) / norm(r)) ^ 2)) * r(2);
                a_j2_z = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (3 - (5 * (r(3) / norm(r)) ^ 2)) * r(3);
                a_j2 = [a_j2_x; a_j2_y; a_j2_z];
                a = a + a_j2;
            end

            r_dot = v;
            v_dot = a;

            statedot = [r_dot; v_dot];
        end

        function statedot = dynamics_relative_deputy(t, t_history, state_d, state_c_history, simulation_settings)

            % State input is relative state in rtn, chief initial state in 
            rho = state_d(1:3);
            drho = state_d(4:6);

            % Extract Closest State
            [~, idx] = min(abs(t_history - t));
            closest_state = state_c_history(idx,:)';
            r_0 = norm(closest_state(1:3));

            h = norm(cross(closest_state(1:3), closest_state(4:6)));
            omega = h / r_0^2;
            omega_dot = -2 * h * (dot(closest_state(1:3), closest_state(4:6)) / norm(closest_state(1:3))) / (r_0 ^ 3);
        
            k = -constants.mu / ((r_0 + rho(1))^2 + rho(2)^2 + rho(3)^2)^(3/2);
            xddot = k * (r_0 + rho(1)) + (constants.mu / r_0^2) + (2 * omega * drho(2)) + (omega_dot * rho(2)) + (rho(1) * omega^2);
            yddot = (k * rho(2)) - (2 * omega * drho(1)) - (omega_dot * rho(1)) + (rho(2) * omega^2);
            zddot = k * rho(3);

            ddrho = [xddot; yddot; zddot];

            statedot = [drho; ddrho];
        end

        function [state_history_kep, t_kep] = propagate_keplerian_orbit(a, e, i, RAAN, w, v, time_span, dt)
            % Propagate orbit using Keplerian motion
            n = sqrt(constants.mu / a ^ 3);
            state_history_kep = zeros(length(time_span), 6);
            
            for j = 1:length(time_span)
                state_kep = util.OE2ECI([a, e, i, RAAN, w, v]);
                state_history_kep(j,:) = state_kep';
                
                % Update true anomaly for next time step
                E = 2 * atan2(tan(v/2) * sqrt((1 - e) / (1 + e)), 1);
                M = E - e*sin(E);
                M_new = M + (n * dt);
                E_new = util.MtoE(M_new, e, 10^(-9));
                v_new = 2 * atan2(sqrt((1 + e) / (1 - e)) * tan(E_new / 2), 1);
                v = v_new;
            end
            
            t_kep = time_span;
        end

        function state_rtn = HCW_propogation(t, initial_conditions_chief_oes, initial_conditions_deputy_rtn)
            a = initial_conditions_chief_oes(1);
            n = sqrt(constants.mu / (a^3));
            a_matrix = [a, 0, 0, 0, 0, 0;
                        0, a, 0, 0, 0, 0;
                        0, 0, a, 0, 0, 0;
                        0, 0, 0, a*n, 0, 0;
                        0, 0, 0, 0, a*n, 0;
                        0, 0, 0, 0, 0, a*n];

            integration_constants = inv(util.calculate_hcw_matrix(0, initial_conditions_chief_oes)) * inv(a_matrix) * initial_conditions_deputy_rtn;
            state_rtn = a_matrix * util.calculate_hcw_matrix(t, initial_conditions_chief_oes) * integration_constants;
        end
    end
end