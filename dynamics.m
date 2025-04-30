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

        function state_rtn = YA_propogation(f, t, initial_conditions_chief_oes, initial_conditions_deputy_rtn)
            a = initial_conditions_chief_oes(1);
            e = initial_conditions_chief_oes(2);
            n = sqrt(constants.mu / (a^3));
            nu = sqrt(1 - e^2);
            a_matrix = [a * nu^2, 0, 0, 0, 0, 0;
                        0, a * nu^2, 0, 0, 0, 0;
                        0, 0, a * nu^2, 0, 0, 0;
                        0, 0, 0, a*n/nu, 0, 0;
                        0, 0, 0, 0, a*n/nu, 0;
                        0, 0, 0, 0, 0, a*n/nu];
            integration_constants = inv(util.calculate_ya_matrix(initial_conditions_chief_oes(6), 0, initial_conditions_chief_oes)) * inv(a_matrix) * initial_conditions_deputy_rtn;
            state_rtn = a_matrix * util.calculate_ya_matrix(f, t, initial_conditions_chief_oes) * integration_constants;
        end

        function state_rtn = propagate_with_roe_circular(t, initial_conditions_chief_oes, initial_roe, current_chief_oe)
            a = initial_conditions_chief_oes(1);
            e = initial_conditions_chief_oes(2);

            da = initial_roe(1);
            dl = initial_roe(2);
            de_x = initial_roe(3);
            de_y = initial_roe(4);
            di_x = initial_roe(5);
            di_y = initial_roe(6);

            de = sqrt(de_x^2 + de_y^2);
            di = sqrt(di_x^2 + di_y^2);
            
            phi_e = atan2(de_y, de_x);
            theta_i = atan2(di_y, di_x);
        
            u = current_chief_oe(5) + util.TtoM(current_chief_oe(6), e);
            u0 = initial_conditions_chief_oes(5) + util.TtoM(initial_conditions_chief_oes(6), e);
            n = sqrt(constants.mu/a^3);
            
            delta_u = u - u0 + 2*pi*floor(n*t/(2*pi));

            dr_r = (da - de * cos(u - phi_e));
            dr_t = (dl - 1.5 * da * delta_u + 2 * de * sin(u - phi_e));
            dr_n = (di * sin(u - theta_i));

            dv_r = n * (de * sin(u - phi_e));
            dv_t = n * (-1.5 * da + 2 * de * cos(u - phi_e));
            dv_n = n * (di * cos(u - theta_i));
        
            state_rtn = [dr_r; dr_t; dr_n; dv_r; dv_t; dv_n];
        end

        function state_rtn = propagate_with_roe_eccentric(t, initial_conditions_chief_oes, initial_roe, current_chief_oe)
            a = initial_conditions_chief_oes(1);
            e = initial_conditions_chief_oes(2);
            i = initial_conditions_chief_oes(3);
            
            da = initial_roe(1);
            dl0 = initial_roe(2);
            de_x = initial_roe(3);
            de_y = initial_roe(4);
            di_x = initial_roe(5);
            di_y = initial_roe(6);
            
            f = current_chief_oe(6);
            w = current_chief_oe(5);
            u = w + f;
            
            e_x = e * cos(w);
            e_y = e * sin(w);
            
            n = sqrt(constants.mu/a^3);
            eta = sqrt(1 - e^2);
            k = 1 + e*cos(f);
            k_prime = -e*sin(f);
            
            n_t = n * t;
            
            b_x_1 = (1/k) + (3/2)*k_prime*(n_t/eta^3);
            b_x_2 = -k_prime/eta^3;
            b_x_3 = (1/eta^3) * (e_x*((k-1)/(1+eta)) - cos(u));
            b_x_4 = (1/eta^3) * (e_y*((k-1)/(1+eta)) - sin(u));
            b_x_6 = (k_prime/eta^3) * cot(i);
            
            b_y_1 = -(3/2)*k*(n_t/eta^3);
            b_y_2 = k/eta^3;
            b_y_3 = (1/eta^2) * ((1 + (1/k))*sin(u) + (e_y/k) + (k/eta)*(e_y/(1+eta)));
            b_y_4 = -(1/eta^2) * ((1 + (1/k))*cos(u) + (e_x/k) + (k/eta)*(e_x/(1+eta)));
            b_y_6 = ((1/k) - (k/eta^3)) * cot(i);
            
            b_z_5 = (1/k) * sin(u);
            b_z_6 = -(1/k) * cos(u);
            
            b_xdot_1 = (k_prime/2) + (3/2)*k^2*(1-k)*(n_t/eta^3);
            b_xdot_2 = -k^2*(k-1)/eta^3;
            b_xdot_3 = (k^2/eta^3) * (eta*sin(u) + e_y*((k-1)/(1+eta)));
            b_xdot_4 = -(k^2/eta^3) * (eta*cos(u) - e_x*((k-1)/(1+eta)));
            b_xdot_6 = -(k^2/eta^3) * (k-1) * cot(i);
            
            b_ydot_1 = -(3/2)*k * (1 + k_prime*(n_t/eta^3));
            b_ydot_2 = (k^2/eta^3) * k_prime;
            b_ydot_3 = (1 + (k^2/eta^3))*cos(u) + e_x*(k/eta^2)*(1 + (k/eta)*((1-k)/(1+eta)));
            b_ydot_4 = (1 + (k^2/eta^3))*sin(u) + e_y*(k/eta^2)*(1 + (k/eta)*((1-k)/(1+eta)));
            b_ydot_6 = -(1 + (k^2/eta^3)) * k_prime * cot(i);
            
            b_zdot_5 = cos(u) + e_x;
            b_zdot_6 = sin(u) + e_y;
            
            dr_r = a * eta^2 * (b_x_1 * da + b_x_2 * dl0 + b_x_3 * de_x + b_x_4 * de_y + 0 * di_x + b_x_6 * di_y);
            dr_t = a * eta^2 * (b_y_1 * da + b_y_2 * dl0 + b_y_3 * de_x + b_y_4 * de_y + 0 * di_x + b_y_6 * di_y);
            dr_n = a * eta^2 * (0 * da + 0 * dl0 + 0 * de_x + 0 * de_y + b_z_5 * di_x + b_z_6 * di_y);
            
            dv_r = a * n / eta * (b_xdot_1 * da + b_xdot_2 * dl0 + b_xdot_3 * de_x + b_xdot_4 * de_y + 0 * di_x + b_xdot_6 * di_y);
            dv_t = a * n / eta * (b_ydot_1 * da + b_ydot_2 * dl0 + b_ydot_3 * de_x + b_ydot_4 * de_y + 0 * di_x + b_ydot_6 * di_y);
            dv_n = a * n / eta * (0 * da + 0 * dl0 + 0 * de_x + 0 * de_y + b_zdot_5 * di_x + b_zdot_6 * di_y);
            
            state_rtn = [dr_r; dr_t; dr_n; dv_r; dv_t; dv_n];
        end
    end
end