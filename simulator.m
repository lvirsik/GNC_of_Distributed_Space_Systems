classdef simulator
    properties
        initial_conditions
        initial_conditions_deputy
        time_span
        dt
        simulation_settings
        graphics_settings
    end

    methods
        function self = simulator(initial_conditions_chief, initial_conditions_deputy, num_orbits, time_step, simulation_settings, graphics_settings)
            self.initial_conditions = initial_conditions_chief;
            self.initial_conditions_deputy = initial_conditions_deputy;
            time = (num_orbits * (2*pi*sqrt(self.initial_conditions(1)^3 / constants.mu)));
            time_span = 0:time_step:time;
            self.dt = time_step;
            self.time_span = time_span;
            self.simulation_settings = simulation_settings;
            self.graphics_settings = graphics_settings;
        end

        function run_simulator(self)
            a = self.initial_conditions(1);
            e = self.initial_conditions(2);
            i = deg2rad(self.initial_conditions(3));
            RAAN = deg2rad(self.initial_conditions(4));
            w = deg2rad(self.initial_conditions(5));
            v = deg2rad(self.initial_conditions(6));
            result.initial_conditions = self.initial_conditions;
            result.dt = self.dt;

            if self.simulation_settings.numerical_propogation
                initial_state_eci = util.OE2ECI(a, e, i, RAAN, w, v);

                % Perform Simulation
                options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
                if self.simulation_settings.relative_deputy
                    initial_state = [self.initial_conditions_deputy; initial_state_eci];
                    [t, state_history] = ode45(@(t, state_history) dynamics.dynamics_with_relative(t, state_history, self.simulation_settings), self.time_span, initial_state, options);
                    result.state_history_num = state_history(:, 7:12);
                    result.relative_state_history = state_history(:, 1:6);
                    result.t_num = t;
                else
                    [t, state_history_num] = ode45(@(t, state_history_num) dynamics.two_body_dynamics(t, state_history_num, self.simulation_settings), self.time_span, initial_state_eci, options);
                    result.state_history_num = state_history_num;
                    result.t_num = t;
                end

                % Implementation for part (c) - compute deputy orbit from chief orbit
                if self.simulation_settings.absolute_deputy
                    % Create initial state for deputy by applying variations to chief
                    % Convert relative state in RTN to ECI
                    rho_RTN = self.initial_conditions_deputy(1:3);
                    drho_RTN = self.initial_conditions_deputy(4:6);
                    
                    % Calculate RTN basis vectors directly
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
                    
                    % Propagate deputy orbit using the existing two_body_dynamics
                    [t_deputy, state_history_deputy] = ode45(@(t, state) dynamics.two_body_dynamics(t, state, self.simulation_settings), self.time_span, deputy_initial_state_eci, options);
                    
                    % Transform deputy orbit to RTN frame relative to chief
                    deputy_in_rtn = zeros(length(t), 6);
                    for j = 1:length(t)
                        % Get chief state at this time point
                        r_chief_j = result.state_history_num(j, 1:3)';
                        v_chief_j = result.state_history_num(j, 4:6)';
                        
                        % Calculate RTN basis vectors
                        R_hat_j = r_chief_j / norm(r_chief_j);
                        h_vec_j = cross(r_chief_j, v_chief_j);
                        N_hat_j = h_vec_j / norm(h_vec_j);
                        T_hat_j = cross(N_hat_j, R_hat_j);
                        
                        % ECI to RTN transformation matrix
                        R_eci2rtn_j = [R_hat_j, T_hat_j, N_hat_j]';
                        
                        % Deputy position and velocity in ECI
                        r_deputy_j = state_history_deputy(j, 1:3)';
                        v_deputy_j = state_history_deputy(j, 4:6)';
                        
                        % Calculate relative position in RTN
                        rho_j = R_eci2rtn_j * (r_deputy_j - r_chief_j);
                        
                        % Calculate relative velocity in RTN
                        omega_j = norm(h_vec_j) / (norm(r_chief_j)^2);
                        omega_vec_j = omega_j * N_hat_j;
                        drho_j = R_eci2rtn_j * (v_deputy_j - v_chief_j - cross(omega_vec_j, r_deputy_j - r_chief_j));
                        
                        deputy_in_rtn(j, :) = [rho_j; drho_j]';
                    end
                    
                    result.deputy_state_history = state_history_deputy;
                    result.deputy_in_rtn = deputy_in_rtn;
                    result.t_deputy = t_deputy;
                end

                % Gather orbital element history and other interesting information
                oe_history = zeros(length(result.state_history_num), 6);
                ecc_vector_history = zeros(length(result.state_history_num), 3);
                ang_mom_history = zeros(length(result.state_history_num), 3);
                energy_history = zeros(length(result.state_history_num), 1);
                for j = 1:length(result.state_history_num)
                    oe_history(j, :) = util.ECI2OE(result.state_history_num(j, :));
                    ecc_vector_history(j, :) = util.get_ecc_vector(result.state_history_num(j, :));
                    ang_mom_history(j, :) = util.get_ang_momentum(result.state_history_num(j, :));
                    energy_history(j, :) = util.get_energy(result.state_history_num(j, :));
                end
                result.oe_history = oe_history;
                result.ecc_vector_history = ecc_vector_history;
                result.ang_mom_history = ang_mom_history;
                result.energy_history = energy_history;
            end

            if self.simulation_settings.keplerian_propogation
                n = sqrt(constants.mu / a ^ 3);
                state_history_kep = zeros(length(self.time_span), 6);
                for j=1:length(self.time_span)
                    state_kep = util.OE2ECI(a, e, i, RAAN, w, v);
                    state_history_kep(j,:) = state_kep';
                    E = 2 * atan2(tan(v/2) * sqrt((1 - e) / (1 + e)), 1);
                    M = E - e*sin(E);
                    M_new = M + (n * self.dt);
                    E_new = util.MtoE(M_new, e, 10^(-9));
                    v_new = 2 * atan2(sqrt((1 + e) / (1 - e)) * tan(E_new / 2), 1);
                    v = v_new;
                end
                result.state_history_kep = state_history_kep;
                result.t_kep = self.time_span;
            end
            plotter(result, self.graphics_settings);
        end
    end
end