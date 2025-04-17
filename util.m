classdef util
    methods (Static)
        function state_eci = OE2ECI(a, e, i, RAAN, w, v)
            p = a * (1 - e^2);
            r = p / (1 + e*cos(v));
            rPQW = [r*cos(v); r*sin(v); 0];
            vPQW = [sqrt(constants.mu/p) * -sin(v); sqrt(constants.mu/p)*(e + cos(v)); 0];
        
            R1 = [cos(-RAAN) sin(-RAAN) 0;...
                -sin(-RAAN) cos(-RAAN) 0;...
                0       0       1];
            R2 = [1  0        0;...
                0  cos(-i) sin(-i);...
                0 -sin(-i) cos(-i)];
            R3 = [cos(-w) sin(-w) 0;...
                -sin(-w) cos(-w) 0;...
                0        0        1];
            R = R1 * R2 * R3;
        
            rECI = R * rPQW;
            vECI = R * vPQW;
            state_eci = [rECI; vECI];
        end

        function H = get_ang_momentum(state)
            R = state(1:3);
            V = state(4:6);
            H = cross(R, V);
        end

        function E = get_ecc_vector(state)
            R = state(1:3);
            V = state(4:6);
            H = util.get_ang_momentum(state);
            mu = constants.mu;
            E = (cross(V,H)/mu)-(R/norm(R));
        end

        function E = get_energy(state)
            oes = util.ECI2OE(state);
            a = oes(1);
            E = -constants.mu / (2 * a);
        end

        function OEs = ECI2OE(state)
            mu = constants.mu;
            R = state(1:3);
            V = state(4:6);
            r = norm(R);
            
            %Get to Perifocal
            H = util.get_ang_momentum(state);
            h = norm(H);
            E = util.get_ecc_vector(state);
            e = norm(E);
            p = (h)^2/mu;
            a = p / (1 - e^2);
            n = sqrt(mu/a^3);
            E = atan2((dot(R, V)/(n*a^2)), (1-(r/a)));
            v = 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2));
            
            %Find Angles
            W = [H(1)/h H(2)/h H(3)/h];
            i = atan2(sqrt(W(1)^2 + W(2)^2), W(3));
            RAAN = atan2(W(1), -W(2));
            w = atan2((R(3)/sin(i)), R(1)*cos(RAAN)+R(2)*sin(RAAN)) - v;
            if rad2deg(w) > 180
                w = w - 2*pi;
            end
            
            OEs = [a, e, i, RAAN, w, v];
        end

        function state_pqw = OE2PQW(a, e, i, RAAN, w, v)
            p = a * (1 - e^2);
            r = p / (1 + e*cos*(v));
            P = r*cos(v);
            Q = r*sin(v);
            W = 0;
            rPQW = [P; Q; W];
            vPQW = sqrt(constants.mu / p) * [-sin(v); e + cos(v); 0];
            state_pqw = [rPQW; vPQW];
        end

        function E = MtoE(M, e, err)
            E = M; N = 0;
            delta = (E - e*sin(E) - M) / (1 - e*cos(E));
            while (delta > err || N < 100)
                delta = (E - e*sin(E) - M) / (1 - e*cos(E));
                E = E - delta;
                N = N + 1;
            end
        end

        function state_rtn = ECI2RTN(state_eci, reference_state_eci)
            r_eci = reference_state_eci(1:3);
            v_eci = reference_state_eci(4:6);
            n = cross(r_eci, v_eci);

            R = r_eci / norm(r_eci);
            N = n / norm(n);
            T = cross(N, R);

            R_eci2rtn = [R, T, N]';
            diff = state_eci(1:3) - reference_state_eci(1:3);

            r_rtn = R_eci2rtn * diff;
            
            f_dot = norm(cross(r_eci, v_eci)) / norm(r_eci)^2;
            w = [0,0,f_dot]';
            v_rtn = R_eci2rtn * (state_eci(4:6) - reference_state_eci(4:6)) - cross(w, r_rtn);
            state_rtn = [r_rtn; v_rtn];
        end

        function R_eci2rtn = R_ECI2RTN(reference_state_eci)
            r_eci = reference_state_eci(1:3)';
            v_eci = reference_state_eci(4:6)';
            n = cross(r_eci, v_eci);

            R = r_eci / norm(r_eci);
            N = n / norm(n);
            T = cross(N, R);

            R_eci2rtn = [R, T, N]';
        end

        function state_eci = RTN2ECI(state_rtn, reference_state_eci)
            rho_RTN = state_rtn(1:3);
            drho_RTN = state_rtn(4:6);
            r_chief = reference_state_eci(1:3);
            v_chief = reference_state_eci(4:6);
            R_hat = r_chief / norm(r_chief);
            h_vec = cross(r_chief, v_chief);
            N_hat = h_vec / norm(h_vec);
            T_hat = cross(N_hat, R_hat);
            R_rtn2eci = [R_hat, T_hat, N_hat];
            
            % Calculate deputy position in ECI
            r_deputy_eci = r_chief + R_rtn2eci * rho_RTN;
            
            % Calculate angular velocity of RTN frame
            omega = norm(h_vec) / (norm(r_chief)^2);
            omega_vec = omega * N_hat;
            
            % Calculate deputy velocity in ECI
            v_deputy_eci = v_chief + R_rtn2eci * drho_RTN + cross(omega_vec, R_rtn2eci * rho_RTN);
            
            % Combine to get deputy initial state in ECI
            state_eci = [r_deputy_eci; v_deputy_eci];
        end

        function state_rtn_history = ECI2RTN_history(deputy_state_history, chief_state_history)
            % Transform deputy orbit to RTN frame relative to chief
            state_rtn_history = zeros(length(deputy_state_history), 6);
            for j = 1:length(state_rtn_history)           
                state_rtn_history(j, :) = util.ECI2RTN(deputy_state_history(j, :)', chief_state_history(j,:)');
            end
        end

        function state_eci_history = RTN2ECI_history(deputy_state_history, chief_state_history)
            % Transform deputy orbit to RTN frame relative to chief
            state_eci_history = zeros(length(deputy_state_history), 6);
            for j = 1:length(state_eci_history)         
                state_eci_history(j, :) = util.RTN2ECI(deputy_state_history(j, :)', chief_state_history(j,:)');
            end
        end
        
        function [oe_history, ecc_vector_history, ang_mom_history, energy_history] = calculate_orbit_history(state_history)
            % Calculate orbital elements history and other orbital parameters
            n_steps = length(state_history);
            oe_history = zeros(n_steps, 6);
            ecc_vector_history = zeros(n_steps, 3);
            ang_mom_history = zeros(n_steps, 3);
            energy_history = zeros(n_steps, 1);
            
            for j = 1:n_steps
                oe_history(j, :) = util.ECI2OE(state_history(j, :));
                ecc_vector_history(j, :) = util.get_ecc_vector(state_history(j, :));
                ang_mom_history(j, :) = util.get_ang_momentum(state_history(j, :));
                energy_history(j, :) = util.get_energy(state_history(j, :));
            end
        end

        function [delta_v, optimal_time_index, maneuver_point] = calculate_drift_correction(deputy_state_history, chief_state_history, t)
            % Calculate orbital elements for both satellites
            [deputy_oe, ~, ~, ~] = util.calculate_orbit_history(deputy_state_history);
            [chief_oe, ~, ~, ~] = util.calculate_orbit_history(chief_state_history);
            
            % Calculate semi-major axis difference
            delta_a = deputy_oe(1, 1) - chief_oe(1, 1);

            % Find points of minimum and maximum radius (approximate apogee/perigee)
            deputy_radius = zeros(length(t), 1);
            for i = 1:length(t)
                deputy_radius(i) = norm(deputy_state_history(i, 1:3));
            end
            [~, min_r_idx] = min(deputy_radius(2:end));
            [~, max_r_idx] = max(deputy_radius(2:end));
            
            % Calculate values at these points
            r_min = deputy_radius(min_r_idx);
            r_max = deputy_radius(max_r_idx);
            v_min = norm(deputy_state_history(min_r_idx, 4:6));
            v_max = norm(deputy_state_history(max_r_idx, 4:6));
            
            % Calculate required velocities for bounded motion
            a_target = chief_oe(1, 1); % Target semi-major axis = chief's
            v_req_min = sqrt(constants.mu * (2/r_min - 1/a_target));
            v_req_max = sqrt(constants.mu * (2/r_max - 1/a_target));
            
            % Calculate delta-v at both locations
            dv_min_r = v_req_min - v_min;
            dv_max_r = v_req_max - v_max;
            
            % Choose the more efficient maneuver
            if abs(dv_min_r) <= abs(dv_max_r)
                delta_v = dv_min_r;
                optimal_time_index = min_r_idx;
                maneuver_point = [r_min, t(min_r_idx)];
            else
                delta_v = dv_max_r;
                optimal_time_index = max_r_idx;
                maneuver_point = [r_max, t(max_r_idx)];
            end
        end
    end
end