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
            v_rtn = R_eci2rtn * state_eci(4:6) - cross(w, r_rtn);
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
        
        function [state_history_kep, t_kep] = propagate_keplerian_orbit(a, e, i, RAAN, w, v, time_span, dt)
            % Propagate orbit using Keplerian motion
            n = sqrt(constants.mu / a ^ 3);
            state_history_kep = zeros(length(time_span), 6);
            
            for j = 1:length(time_span)
                state_kep = util.OE2ECI(a, e, i, RAAN, w, v);
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
        
        function [deputy_initial_state_eci, deputy_in_rtn, deputy_state_history, t_deputy] = propagate_deputy_orbit(initial_conditions_deputy, initial_state_eci, time_span, chief_state_history, t, simulation_settings)
            % Create initial state for deputy by applying variations to chief
            % Convert relative state in RTN to ECI
            rho_RTN = initial_conditions_deputy(1:3);
            drho_RTN = initial_conditions_deputy(4:6);
            
            % Calculate RTN basis vectors
            r_chief = initial_state_eci(1:3);
            v_chief = initial_state_eci(4:6);
            
            R_hat = r_chief / norm(r_chief);
            h_vec = cross(r_chief, v_chief);
            N_hat = h_vec / norm(h_vec);
            T_hat = cross(N_hat, R_hat);
            
            % RTN to ECI transformation matrix (each column is a basis vector)
            R_rtn2eci = [R_hat, T_hat, N_hat];
            
            % Calculate deputy position in ECI
            r_deputy_eci = r_chief + R_rtn2eci * rho_RTN;
            
            % Calculate angular velocity of RTN frame
            omega = norm(h_vec) / (norm(r_chief)^2);
            omega_vec = omega * N_hat;
            
            % Calculate deputy velocity in ECI
            v_deputy_eci = v_chief + R_rtn2eci * drho_RTN + cross(omega_vec, R_rtn2eci * rho_RTN);
            
            % Combine to get deputy initial state in ECI
            deputy_initial_state_eci = [r_deputy_eci; v_deputy_eci];
            
            % Propagate deputy orbit using ODE45
            options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
            [t_deputy, deputy_state_history] = ode45(@(t, state) dynamics.two_body_dynamics(t, state, simulation_settings), time_span, deputy_initial_state_eci, options);
            
            % Transform deputy orbit to RTN frame relative to chief
            deputy_in_rtn = util.transform_deputy_to_rtn(deputy_state_history, chief_state_history, t);
        end
        
        function deputy_in_rtn = transform_deputy_to_rtn(deputy_state_history, chief_state_history, t)
            % Transform deputy orbit to RTN frame relative to chief
            deputy_in_rtn = zeros(length(t), 6);
            
            for j = 1:length(t)
                % Get chief state at this time point
                r_chief_j = chief_state_history(j, 1:3)';
                v_chief_j = chief_state_history(j, 4:6)';
                
                % Calculate RTN basis vectors
                R_hat_j = r_chief_j / norm(r_chief_j);
                h_vec_j = cross(r_chief_j, v_chief_j);
                N_hat_j = h_vec_j / norm(h_vec_j);
                T_hat_j = cross(N_hat_j, R_hat_j);
                
                % ECI to RTN transformation matrix
                R_eci2rtn_j = [R_hat_j, T_hat_j, N_hat_j]';
                
                % Deputy position and velocity in ECI
                r_deputy_j = deputy_state_history(j, 1:3)';
                v_deputy_j = deputy_state_history(j, 4:6)';
                
                % Calculate relative position in RTN
                rho_j = R_eci2rtn_j * (r_deputy_j - r_chief_j);
                
                % Calculate relative velocity in RTN
                omega_j = norm(h_vec_j) / (norm(r_chief_j)^2);
                omega_vec_j = omega_j * N_hat_j;
                drho_j = R_eci2rtn_j * (v_deputy_j - v_chief_j - cross(omega_vec_j, r_deputy_j - r_chief_j));
                
                deputy_in_rtn(j, :) = [rho_j; drho_j]';
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
            deputy_oe = zeros(length(t), 6);
            chief_oe = zeros(length(t), 6);
            
            for i = 1:length(t)
                deputy_oe(i, :) = util.ECI2OE(deputy_state_history(i, :));
                chief_oe(i, :) = util.ECI2OE(chief_state_history(i, :));
            end
            
            % Calculate semi-major axis difference
            delta_a = deputy_oe(:, 1) - chief_oe(:, 1);
            disp(['Current semi-major axis difference: ', num2str(delta_a(end)), ' meters']);
            
            % Find points of minimum and maximum radius (approximate apogee/perigee)
            deputy_radius = zeros(length(t), 1);
            for i = 1:length(t)
                deputy_radius(i) = norm(deputy_state_history(i, 1:3));
            end
            
            [~, min_r_idx] = min(deputy_radius);
            [~, max_r_idx] = max(deputy_radius);
            
            % Calculate velocities at these points
            v_min_r = norm(deputy_state_history(min_r_idx, 4:6));
            v_max_r = norm(deputy_state_history(max_r_idx, 4:6));
            
            % Calculate required velocities for bounded motion
            a_target = chief_oe(1, 1); % Target semi-major axis = chief's
            r_min = deputy_radius(min_r_idx);
            r_max = deputy_radius(max_r_idx);
            
            v_req_min_r = sqrt(constants.mu * (2/r_min - 1/a_target));
            v_req_max_r = sqrt(constants.mu * (2/r_max - 1/a_target));
            
            % Calculate delta-v at both locations
            dv_min_r = abs(v_req_min_r - v_min_r);
            dv_max_r = abs(v_req_max_r - v_max_r);
            
            % Choose the more efficient maneuver
            if dv_min_r <= dv_max_r
                delta_v = dv_min_r;
                optimal_time_index = min_r_idx;
                maneuver_point = [r_min, t(min_r_idx)];
                disp(['Optimal maneuver at minimum radius (perigee), delta-v = ', num2str(delta_v), ' m/s']);
            else
                delta_v = dv_max_r;
                optimal_time_index = max_r_idx;
                maneuver_point = [r_max, t(max_r_idx)];
                disp(['Optimal maneuver at maximum radius (apogee), delta-v = ', num2str(delta_v), ' m/s']);
            end
        end
        
        function [t_combined, deputy_state_combined, deputy_rtn_combined] = apply_maneuver(result, optimal_time_index, time_span, simulation_settings)
            % Get the state at maneuver point
            maneuver_state = result.deputy_state_history(optimal_time_index, :);
            
            % Calculate the new velocity after maneuver
            v_deputy = maneuver_state(4:6);
            v_unit = v_deputy / norm(v_deputy);
            
            % We need to match chief's semi-major axis
            a_chief = util.ECI2OE(result.state_history_num(optimal_time_index, :));
            a_chief = a_chief(1);
            r_deputy = norm(maneuver_state(1:3));
            
            % Calculate required velocity magnitude for the desired orbit
            v_new_mag = sqrt(constants.mu * (2/r_deputy - 1/a_chief));
            v_new = v_unit * v_new_mag;
            
            % Apply the maneuver
            post_maneuver_state = maneuver_state;
            post_maneuver_state(4:6) = v_new;
            
            % Continue propagation from maneuver point
            t_remaining = time_span(time_span > result.t_num(optimal_time_index));
            
            if ~isempty(t_remaining)
                options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
                [t_post, state_post] = ode45(@(t, state) dynamics.two_body_dynamics(t, state, simulation_settings), t_remaining, post_maneuver_state, options);
                
                % Store post-maneuver trajectory
                t_combined = [result.t_num(1:optimal_time_index); t_post];
                deputy_state_combined = [result.deputy_state_history(1:optimal_time_index, :); state_post];
                
                % Transform post-maneuver trajectory to RTN
                deputy_rtn_post = util.transform_post_maneuver_to_rtn(state_post, t_post, result);
                
                % Combine trajectories
                deputy_rtn_combined = [result.deputy_in_rtn(1:optimal_time_index, :); deputy_rtn_post];
            else
                t_combined = result.t_num(1:optimal_time_index);
                deputy_state_combined = result.deputy_state_history(1:optimal_time_index, :);
                deputy_rtn_combined = result.deputy_in_rtn(1:optimal_time_index, :);
            end
        end
        
        function deputy_rtn_post = transform_post_maneuver_to_rtn(state_post, t_post, result)
            % Transform post-maneuver trajectory to RTN
            deputy_rtn_post = zeros(length(t_post), 6);
            
            for j = 1:length(t_post)
                % Find closest chief state time
                [~, t_idx] = min(abs(result.t_num - t_post(j)));
                
                % Chief state
                r_chief_j = result.state_history_num(t_idx, 1:3)';
                v_chief_j = result.state_history_num(t_idx, 4:6)';
                
                % Calculate RTN basis
                R_hat_j = r_chief_j / norm(r_chief_j);
                h_vec_j = cross(r_chief_j, v_chief_j);
                N_hat_j = h_vec_j / norm(h_vec_j);
                T_hat_j = cross(N_hat_j, R_hat_j);
                
                R_eci2rtn_j = [R_hat_j, T_hat_j, N_hat_j]';
                
                % Deputy state
                r_deputy_j = state_post(j, 1:3)';
                v_deputy_j = state_post(j, 4:6)';
                
                % Calculate relative position
                rho_j = R_eci2rtn_j * (r_deputy_j - r_chief_j);
                
                % Calculate relative velocity
                omega_j = norm(h_vec_j) / (norm(r_chief_j)^2);
                omega_vec_j = omega_j * N_hat_j;
                drho_j = R_eci2rtn_j * (v_deputy_j - v_chief_j - cross(omega_vec_j, r_deputy_j - r_chief_j));
                
                deputy_rtn_post(j, :) = [rho_j; drho_j]';
            end
        end
    end
end