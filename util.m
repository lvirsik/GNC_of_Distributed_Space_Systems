classdef util
    methods (Static)
        function state_eci = OE2ECI(state)
            a = state(1);
            e = state(2);
            i = state(3);
            RAAN = state(4);
            w = state(5);
            v = state(6);

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
            if i == 0
                i = i + 0.0000001;
            end
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
            T = cross(N, R) / norm(cross(N, R));

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
    
        function hcw_matrix = calculate_hcw_matrix(t, oes)
            n = sqrt(constants.mu / (oes(1)^3));
            hcw_matrix = [1, sin(n*t), cos(n*t), 0, 0, 0;
                        (-3/2)*n*t, 2*cos(n*t), -2*sin(n*t), 1, 0, 0;
                        0, 0, 0, 0, sin(n*t), cos(n*t);
                        0, cos(n*t), -sin(n*t), 0, 0, 0;
                        (-3/2), -2*sin(n*t), -2*cos(n*t), 0, 0, 0;
                        0, 0, 0, 0, cos(n*t), -sin(n*t)];
        end

        function [RR, RV, VR, VV] = t2PSI(t, a)
            n = sqrt(constants.mu / a^3);
            
            RR = [4 - 3*cos(n*t),       0,              0;
                    6*(sin(n*t) - n*t),   1,              0;
                    0,                    0,              cos(n*t)];
                
            RV = [sin(n*t)/n,           2*(1 - cos(n*t))/n,     0;
                    2*(cos(n*t) - 1)/n,   (4*sin(n*t) - 3*n*t)/n, 0;
                    0,                    0,                     sin(n*t)/n];
                
            VR = [3*n*sin(n*t),         0,              0;
                    6*n*(cos(n*t) - 1),   0,              0;
                    0,                    0,             -n*sin(n*t)];
                
            VV = [cos(n*t),             2*sin(n*t),            0;
                    -2*sin(n*t),          4*cos(n*t) - 3,        0;
                    0,                    0,                    cos(n*t)];
        end
    
        function ya_matrix = calculate_ya_matrix(f, t, oes)
            e = oes(2);
            a = oes(1);
            n = sqrt(constants.mu / (a^3));
            nu = sqrt(1 - e^2);
            k = 1 + (e * cos(f));
            kd = -e * sin(f);
            tau = n * t / (nu^3);

            ya_matrix = [(1/k) + (3*kd*tau/2), sin(f), cos(f), 0, 0, 0;
                        -3*k*tau/2, (1 + (1/k))*cos(f), -(1 + (1/k))*sin(f), 1/k, 0, 0;
                        0, 0, 0, 0, (1/k)*sin(f), (1/k)*cos(f);
                        (kd/2)-(3*(k^2)*(k-1)*tau/2), (k^2)*cos(f), -(k^2)*sin(f), 0, 0, 0;
                        -(3/2)*(k + (kd*tau*k^2)), -(k^2 + 1)*sin(f), -e-((1+k^2) * cos(f)), -kd, 0, 0;
                        0, 0, 0, 0, e + cos(f), -sin(f)];
        end

        function M = TtoM(nu, e)
            E = 2 * atan(sqrt((1-e)/(1+e)) * tan(nu/2));
            if E < 0
                E = E + 2*pi;
            end
            M = E - e * sin(E);
        end

        function roe = ECI2ROE(chief_state, deputy_state)
            chief_oe = util.ECI2OE(chief_state);
            deputy_oe = util.ECI2OE(deputy_state);
        
            a_c = chief_oe(1);
            e_c = chief_oe(2);
            i_c = chief_oe(3);
            RAAN_c = chief_oe(4);
            w_c = chief_oe(5);
            nu_c = chief_oe(6);
            
            a_d = deputy_oe(1);
            e_d = deputy_oe(2);
            i_d = deputy_oe(3);
            RAAN_d = deputy_oe(4);
            w_d = deputy_oe(5);
            nu_d = deputy_oe(6);
        
            M_c = util.TtoM(nu_c, e_c);
            M_d = util.TtoM(nu_d, e_d);
            
            da = (a_d - a_c) / a_c;
            dl = wrapTo2Pi((M_d + w_d) - (M_c + w_c) + (RAAN_d - RAAN_c)*cos(i_c));

            de_x = e_d*cos(w_d) - e_c*cos(w_c);
            de_y = e_d*sin(w_d) - e_c*sin(w_c);
            
            di_x = i_d - i_c;
            di_y = (RAAN_d - RAAN_c)*sin(i_c);
            
            roe = [da; dl; de_x; de_y; di_x; di_y];
        end

        function state_rtn = ROE2RTN(initial_roe, current_chief_oe)
            a = current_chief_oe(1);
            e = current_chief_oe(2);
            i = current_chief_oe(3);
            
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
            
            n_t = 0;
            
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
            
            dr_r = eta^2 * (b_x_1 * da + b_x_2 * dl0 + b_x_3 * de_x + b_x_4 * de_y + 0 * di_x + b_x_6 * di_y);
            dr_t = eta^2 * (b_y_1 * da + b_y_2 * dl0 + b_y_3 * de_x + b_y_4 * de_y + 0 * di_x + b_y_6 * di_y);
            dr_n = eta^2 * (0 * da + 0 * dl0 + 0 * de_x + 0 * de_y + b_z_5 * di_x + b_z_6 * di_y);
            
            dv_r = n / eta * (b_xdot_1 * da + b_xdot_2 * dl0 + b_xdot_3 * de_x + b_xdot_4 * de_y + 0 * di_x + b_xdot_6 * di_y);
            dv_t = n / eta * (b_ydot_1 * da + b_ydot_2 * dl0 + b_ydot_3 * de_x + b_ydot_4 * de_y + 0 * di_x + b_ydot_6 * di_y);
            dv_n = n / eta * (0 * da + 0 * dl0 + 0 * de_x + 0 * de_y + b_zdot_5 * di_x + b_zdot_6 * di_y);
            
            state_rtn = [dr_r; dr_t; dr_n; dv_r; dv_t; dv_n];
        end

        function q_oe = ECI2qOE(state_eci)
            oes = util.ECI2OE(state_eci);
            a = oes(1);
            e = oes(2);
            i = oes(3);
            RAAN = oes(4);
            w = oes(5);
            f = oes(6);

            l = RAAN + w + util.TtoM(f, e);
            ex = e * cos(w + RAAN);
            ey = e * sin(w + RAAN);
            ix = sin(i)*cos(RAAN);
            iy = sin(i)*sin(RAAN);

            q_oe = [a; l; ex; ey; ix; iy];
        end

        function STM = calculate_stm_j2(oes, t)
            tau = t;
            a = oes(1);
            e = oes(2);
            i = oes(3);
            eta = sqrt(1 - e^2);
            J2 = constants.J2;
            R_E = constants.earth_radius;
            mu = constants.mu;
            kappa = (3/4) * J2 * R_E^2 * sqrt(mu) / (a^(7/2) * eta^4);
            E = 1 + eta;
            F = 4 + 3 * eta;
            G = 1 / eta^2;
            P = 3 * cos(i)^2 - 1;
            Q = 5 * cos(i)^2 - 1;
            R = cos(i);
            S = sin(2 * i);
            T = sin(i)^2;
            U = sin(i);
            V = tan(i/2);
            W = cos(i/2)^2;

            % Time-varying angular quantities
            omega_dot = kappa * Q;
            Omega_dot = -2 * kappa * R;

            omega_f = omega_i + omega_dot * tau;
            Omega_f = Omega_i + Omega_dot * tau;

            % Eccentricity vector components
            e_xi = e * cos(omega_i);
            e_yi = e * sin(omega_i);
            e_xf = e * cos(omega_f);
            e_yf = e * sin(omega_f);

            e_xxi = e * cos(omega_i + Omega_i);
            e_yyi = e * sin(omega_i + Omega_i);
            e_xxf = e * cos(omega_f + Omega_f);
            e_yyf = e * sin(omega_f + Omega_f);

            % Inclination vector components
            t_xi = tan(i/2) * cos(Omega_i);
            t_yi = tan(i/2) * sin(Omega_i);
            t_xf = tan(i/2) * cos(Omega_f);
            t_yf = tan(i/2) * sin(Omega_f);

            % Initialize 7x7 identity STM matrix
            Phi = sym(eye(7));

            % Row 2
            Phi(2,1) = -(3/2)*n - (7/2)*kappa*E*P;
            Phi(2,3) = kappa * e_xi * F * G * P * tau;
            Phi(2,4) = kappa * e_yi * F * G * P * tau;
            Phi(2,6) = -kappa * F * S * tau;

            % Row 3
            Phi(3,1) = (7/2) * kappa * e_yf * Q * tau;
            Phi(3,3) = cos(omega_dot * tau) - 4 * kappa * e_xxi * e_yyf * G * Q * tau;
            Phi(3,4) = sin(omega_dot * tau) - 4 * kappa * e_yyi * e_yyf * G * Q * tau;
            Phi(3,5) = 5 * kappa * e_yf * S * tau;

            % Row 4
            Phi(4,1) = -(7/2) * kappa * e_xf * Q * tau;
            Phi(4,3) = sin(omega_dot * tau) + 4 * kappa * e_yi * e_yf * G * Q * tau;
            Phi(4,4) = cos(omega_dot * tau) + 4 * kappa * e_yi * e_xf * G * Q * tau;
            Phi(4,5) = -5 * kappa * e_xf * S * tau;

            % Row 6
            Phi(6,1) = (7/2) * kappa * S * tau;
            Phi(6,3) = -4 * kappa * e_xi * G * S * tau;
            Phi(6,4) = -4 * kappa * e_yi * G * S * tau;
            Phi(6,5) = 2 * kappa * T * tau;

        end
    end
end