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

                % Implementation for part (2c) - compute deputy orbit from chief orbit
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
                    
                    % Part (2e) and (2f) - Calculate and apply maneuver if requested
                    if isfield(self.simulation_settings, 'calculate_maneuver') && self.simulation_settings.calculate_maneuver
                        % Calculate the optimal maneuver to correct the drift
                        [delta_v, optimal_time_index, maneuver_point] = self.calculate_drift_correction(result.deputy_state_history, result.state_history_num, result.t_num);
                        
                        % Save maneuver information
                        result.maneuver_time = result.t_num(optimal_time_index);
                        result.maneuver_delta_v = delta_v;
                        result.maneuver_time_index = optimal_time_index;
                        result.maneuver_point = maneuver_point;
                        
                        % Apply maneuver if requested
                        if isfield(self.simulation_settings, 'apply_maneuver') && self.simulation_settings.apply_maneuver
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
                            t_remaining = self.time_span(self.time_span > result.t_num(optimal_time_index));
                            if ~isempty(t_remaining)
                                [t_post, state_post] = ode45(@(t, state) dynamics.two_body_dynamics(t, state, self.simulation_settings), t_remaining, post_maneuver_state, options);
                                
                                % Store post-maneuver trajectory
                                t_combined = [result.t_num(1:optimal_time_index); t_post];
                                deputy_combined = [result.deputy_state_history(1:optimal_time_index, :); state_post];
                                
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
                                
                                % Combine trajectories
                                deputy_rtn_combined = [result.deputy_in_rtn(1:optimal_time_index, :); deputy_rtn_post];
                                
                                % Save combined results
                                result.t_combined = t_combined;
                                result.deputy_state_combined = deputy_combined;
                                result.deputy_in_rtn_combined = deputy_rtn_combined;
                            end
                        end
                    end
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
        
        function [delta_v, optimal_time_index, maneuver_point] = calculate_drift_correction(self, deputy_state_history, chief_state_history, t)
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
    end
end