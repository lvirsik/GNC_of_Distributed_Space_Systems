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

        function statedot = dynamics_with_relative(t, state, simulation_settings)
            % State input is relative state in rtn, followed by chief state in eci
            rho = state(1:3);
            drho = state(4:6);
            chief_state_rtn = util.ECI2RTN(state(7:12), state(7:12));
            r_0 = [norm(state(7:9)); 0; 0];
            v_0 = chief_state_rtn(4:6);

            h = norm(cross(r_0, v_0));
            omega = h / norm(r_0)^2;
            omega_dot = -2 * h * v_0(1) / (norm(r_0) ^ 3);
        
            k = -constants.mu / ((norm(r_0) + rho(1))^2 + rho(2)^2 + rho(3)^2)^(3/2);
            xddot = k * (norm(r_0) + rho(1)) + (constants.mu / norm(r_0)^2) + (2 * omega * drho(2)) + (omega_dot * rho(2)) + (rho(1) * omega^2);
            yddot = (k * rho(2)) - (2 * omega * drho(1)) - (omega_dot * rho(1)) + (rho(2) * omega^2);
            zddot = k * rho(3);

            ddrho = [xddot; yddot; zddot];

            chief_dot = dynamics.two_body_dynamics(t, state(7:12), simulation_settings);
            
            statedot = [drho; ddrho; chief_dot];
        end
    end
end