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

        function state_rtn = ECI2RTN(state_eci)
            r_eci = state_eci(1:3)';
            v_eci = state_eci(4:6)';
            n = cross(r_eci, v_eci);

            R = r_eci / norm(r_eci);
            N = n / norm(n);
            T = cross(N, R);

            R_eci2rtn = [R, T, N];

            r_rtn = R_eci2rtn * r_eci;

            f_dot = norm(cross(r_eci, v_eci)) / norm(r_eci)^2;
            w = [0,0,f_dot]';
            v_rtn = (R_eci2rtn * v_eci) - cross(w, r_rtn);

            state_rtn = [r_rtn; v_rtn];
        end
    end
 end