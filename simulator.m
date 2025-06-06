classdef simulator < handle
    properties
        initial_conditions_chief
        initial_conditions_deputy
        time_span
        dt
        simulation_settings
        graphics_settings
        current_step
        previous_time
        P
        estimated_state_history
        covariance_history
        delta_v_tracker
        last_control_time
        acceleration
    end

    methods
        function self = simulator(initial_conditions_chief, initial_conditions_deputy, num_orbits, time_step, simulation_settings, graphics_settings)
            self.initial_conditions_chief = initial_conditions_chief;
            self.initial_conditions_deputy = initial_conditions_deputy;
            time = (num_orbits * (2*pi*sqrt(self.initial_conditions_chief(1)^3 / constants.mu)));
            time_span = 0:time_step:time;
            self.dt = time_step;
            self.time_span = time_span;
            self.simulation_settings = simulation_settings;
            self.graphics_settings = graphics_settings;
            self.current_step = 1;
            self.previous_time = 0;
            self.P = diag([1e4, 1e4, 1e4, 1, 1, 1]);
            self.delta_v_tracker = 0;
            self.last_control_time = -inf;
            self.acceleration = [0;0;0];
        end

        function run_simulator(self)
            % Convert Initial State to Orbital Elements
            a = self.initial_conditions_chief(1);
            e = self.initial_conditions_chief(2);
            incl = deg2rad(self.initial_conditions_chief(3));
            RAAN = deg2rad(self.initial_conditions_chief(4));
            w = deg2rad(self.initial_conditions_chief(5));
            v = deg2rad(self.initial_conditions_chief(6));

            % Save initial conditions to result obj
            result.initial_conditions_chief = self.initial_conditions_chief;
            result.dt = self.dt;

            % Initialize Information
            chief_initial_state_eci = util.OE2ECI([a, e, incl, RAAN, w, v]);
            options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

            % Run Propogation for chief satellite
            [result.t_num, result.chief_history_num] = ode45(@(t, state_history_num) dynamics.two_body_dynamics(t, state_history_num, self.simulation_settings), self.time_span, chief_initial_state_eci, options);            

            if self.simulation_settings.simulation_with_flight_computer
                result.dv = zeros(1,3);
                self.delta_v_tracker = 0;
                self.estimated_state_history = zeros(size(result.chief_history_num));
                self.covariance_history = zeros(length(result.chief_history_num), 6, 6);
                initial_conditions_deputy_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
                options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', self.dt);
                [time_output, deputy_state_history] = ode45(@(t, state) self.wrapper_state_to_stateDot(t, state, result), self.time_span, initial_conditions_deputy_eci, options);

                % Compute velocity magnitude at each time step
                velocity_magnitude = vecnorm(deputy_state_history(:, 4:6), 2, 2);  % 2-norm across rows

                % Plot velocity magnitude vs time
                figure;
                plot(time_output, velocity_magnitude, 'LineWidth', 2);
                xlabel('Time [s]');
                ylabel('‖Velocity‖ [km/s]');
                title('Deputy Velocity Magnitude Over Time');
                grid on;
                
                % Transform deputy orbit to RTN frame relative to chief
                result.absolute_state_history = util.ECI2RTN_history(deputy_state_history, result.chief_history_num);
                result.deputy_state_history_eci = deputy_state_history;  
                
                result.estimated_state_history = util.ECI2RTN_history(self.estimated_state_history, result.chief_history_num);
                result.covariance_history = self.covariance_history;
                result.estimated_state_history(end,:) = result.estimated_state_history(end-1,:);

            elseif self.simulation_settings.manuver_instant
                % Initialize Information
                dv_tracker = 0;
                result.dv = [];
                chief_initial_state_eci = util.OE2ECI([a, e, incl, RAAN, w, v]);
                deputy_initial_state_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
                initial_state_eci = [chief_initial_state_eci; deputy_initial_state_eci];
                options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

                orbit_period = (2*pi*sqrt(a^3 / constants.mu));
                num_orbits = round(self.time_span(end) / orbit_period);

                [new_t, new_history] = ode45(@(t, state_history_num) dynamics.wrapper_two_body_relative(t, state_history_num, self.simulation_settings), self.time_span / 2, initial_state_eci, options);

                result.combined_history = new_history;
                result.t_num = new_t;

                initial_state_eci = result.combined_history(end, :);
                chief_eci = initial_state_eci(1:6)';
                deputy_eci = initial_state_eci(7:12)';
                
                    
                % MANUVER WORK DONE HERE
                deputy_rtn = util.ECI2RTN(deputy_eci, chief_eci);
                ti = self.time_span(end) * 3/4;
                tf = self.time_span(end);

                [rr, rv, vr, vv] = util.t2PSI(tf-ti, a);
                r0 = deputy_rtn(1:3);
                v0 = deputy_rtn(4:6);
                new_v = -pinv(rv) * rr * r0
                deputy_rtn(4:6) = new_v;
                %%%%%%
 
                deputy_eci = util.RTN2ECI(deputy_rtn, chief_eci);
                initial_state_eci(7:12) = deputy_eci;

                [new_t, new_history] = ode45(@(t, state_history_num) dynamics.wrapper_two_body_relative(t, state_history_num, self.simulation_settings), self.time_span / 2, initial_state_eci, options);
                result.combined_history = [result.combined_history; new_history];
                result.t_num = [result.t_num; new_t + (ones(size(new_t)) * result.t_num(end))];
                
                % Post Processing
                result.chief_history_num = result.combined_history(:, 1:6);
                result.deputy_history_num = result.combined_history(:, 7:12);
                result.absolute_state_history = util.ECI2RTN_history(result.deputy_history_num, result.chief_history_num);
                result.deputy_state_history_eci = result.deputy_history_num;
            else
                % Initialize Information
                chief_initial_state_eci = util.OE2ECI([a, e, incl, RAAN, w, v]);
                options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

                % Run Propogation for chief satellite
                [result.t_num, result.chief_history_num] = ode45(@(t, state_history_num) dynamics.two_body_dynamics(t, state_history_num, self.simulation_settings), self.time_span, chief_initial_state_eci, options);

                % If there is a relative deputy, propogate using the relative equations of motion
                if self.simulation_settings.relative_deputy
                    [~, result.relative_state_history] = ode45(@(t, state) dynamics.dynamics_relative_deputy(t, result.t_num, state, result.chief_history_num, self.simulation_settings), self.time_span, self.initial_conditions_deputy, options);
                    result.deputy_state_history_eci = util.RTN2ECI_history(result.relative_state_history, result.chief_history_num);
                end
                
                % If there is an absolute deputy, propogate using the 2 body equations of motion
                if self.simulation_settings.absolute_deputy
                    initial_conditions_deputy_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
                    [~, deputy_state_history] = ode45(@(t, state) dynamics.two_body_dynamics(t, state, self.simulation_settings), self.time_span, initial_conditions_deputy_eci, options);
                    
                    % Transform deputy orbit to RTN frame relative to chief
                    result.absolute_state_history = util.ECI2RTN_history(deputy_state_history, result.chief_history_num);
                    result.deputy_state_history_eci = deputy_state_history;
                end

                % If there is a hcw deputy, propogate using hcw equations
                if self.simulation_settings.hcw_deputy
                    hcw_deputy_state_history = zeros(length(self.time_span), 6);
                    for i = 1:length(self.time_span)
                        hcw_deputy_state_history(i,:) = dynamics.HCW_propogation(self.time_span(i), self.initial_conditions_chief, self.initial_conditions_deputy);
                    end
                    result.hcw_state_history = hcw_deputy_state_history;
                end

                % If there is a ya deputy, propogate using ya equations
                if self.simulation_settings.ya_deputy
                    ya_deputy_state_history = zeros(length(self.time_span), 6);
                    for i = 1:length(self.time_span)
                        chief_oes = util.ECI2OE(result.chief_history_num(i, :));
                        ya_deputy_state_history(i,:) = dynamics.YA_propogation(chief_oes(6), self.time_span(i), self.initial_conditions_chief, self.initial_conditions_deputy);
                    end
                    result.ya_state_history = ya_deputy_state_history;
                end

                % If there is a roe deputy, propagate using ROE
                if self.simulation_settings.roe_circular_deputy
                    roe_circular_deputy_state_history = zeros(length(self.time_span), 6);
                    deputy_initial_state_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
                    initial_roe = util.ECI2ROE(chief_initial_state_eci, deputy_initial_state_eci);
                    for idx = 1:length(self.time_span)
                        chief_oes = util.ECI2OE(result.chief_history_num(idx, :));
                        roe_circular_deputy_state_history(idx,:) = dynamics.propagate_with_roe_circular(self.time_span(idx), [a, e, incl, RAAN, w, v], initial_roe, chief_oes);
                    end
                    result.roe_circular_state_history = roe_circular_deputy_state_history;
                end

                if self.simulation_settings.roe_eccentric_deputy
                    roe_eccentric_deputy_state_history = zeros(length(self.time_span), 6);
                    deputy_initial_state_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
                    initial_roe = util.ECI2ROE(chief_initial_state_eci, deputy_initial_state_eci)
                    for idx = 1:length(self.time_span)
                        chief_oes = util.ECI2OE(result.chief_history_num(idx, :));
                        roe_eccentric_deputy_state_history(idx,:) = dynamics.propagate_with_roe_eccentric(self.time_span(idx), [a, e, incl, RAAN, w, v], initial_roe, chief_oes);
                    end
                    result.roe_eccentric_state_history = roe_eccentric_deputy_state_history;
                end

                % If there is a desire to manuver deputy to have bounded motion, calculate the manuver
                if self.simulation_settings.create_bounded_motion
                    [delta_v, optimal_time_index, maneuver_point] = util.calculate_drift_correction(result.deputy_state_history_eci, result.chief_history_num, result.t_num);
                    
                    result.maneuver_time = result.t_num(optimal_time_index);
                    result.maneuver_delta_v = delta_v;
                    result.maneuver_time_index = optimal_time_index;
                    result.maneuver_point = maneuver_point;

                    t_span_1 = self.time_span(1:optimal_time_index);
                    t_span_2 = self.time_span(optimal_time_index+1:end);   

                    [t1, piece1] = ode45(@(t, state) dynamics.dynamics_relative_deputy(t, result.t_num, state, result.chief_history_num, self.simulation_settings), t_span_1, self.initial_conditions_deputy, options);
                    impulse = piece1(end,:) + [0,0,0,0,delta_v,0];
                    [~, piece2] = ode45(@(t, state) dynamics.dynamics_relative_deputy(t + t1(end), result.t_num, state, result.chief_history_num, self.simulation_settings), t_span_2, impulse, options);

                    result.deputy_manuvered = [piece1; piece2];
                end

                [result.oe_history, result.ecc_vector_history, result.ang_mom_history, result.energy_history] = util.calculate_orbit_history(result.chief_history_num);
            end

            % Run keplarian propogator 
            if self.simulation_settings.keplerian_propogation
                [result.state_history_kep, result.t_kep] = dynamics.propagate_keplerian_orbit(a, e, i, RAAN, w, v, self.time_span, self.dt);
            end
            
            % Run Plotter
            plotter(result, self.graphics_settings);
        end

        function continous_control(self)
            % Initialize Information
            dv_tracker = 0;
            result.dv = [];
            a = self.initial_conditions_chief(1);
            chief_initial_state_eci = util.OE2ECI(self.initial_conditions_chief);
            deputy_initial_state_eci = util.RTN2ECI(self.initial_conditions_deputy, chief_initial_state_eci);
            initial_state_eci = [chief_initial_state_eci; deputy_initial_state_eci];
            options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

            orbit_period = (2*pi*sqrt(a^3 / constants.mu));
            num_orbits = round(self.time_span(end) / orbit_period);
            manuvers_per_orbit = 25;
            for i=1:num_orbits*manuvers_per_orbit
                [new_t, new_history] = ode45(@(t, state_history_num) dynamics.wrapper_two_body_relative(t, state_history_num, self.simulation_settings), self.time_span / (num_orbits*manuvers_per_orbit), initial_state_eci, options);

                if i == 1
                    result.combined_history = new_history;
                    result.t_num = new_t;
                else
                    result.combined_history = [result.combined_history; new_history];
                    result.t_num = [result.t_num; new_t + (ones(size(new_t)) * result.t_num(end))];
                end

                initial_state_eci = result.combined_history(end, :);

                chief_eci = initial_state_eci(1:6)';
                deputy_eci = initial_state_eci(7:12)';
                deputy_rtn = util.ECI2RTN(deputy_eci, chief_eci);
                
                result.dv = [result.dv; dv_tracker*ones(size(new_t))];
                burn_dv = 0.01;
                old_v = deputy_rtn(4:6);
                deputy_rtn(4:6) = burn_dv * -deputy_rtn(1:3) / norm(deputy_rtn(1:3));
                dv_tracker = dv_tracker + norm(deputy_rtn(4:6) - old_v);

                deputy_eci = util.RTN2ECI(deputy_rtn, chief_eci);
                initial_state_eci(7:12) = deputy_eci;
            end
            
            % Post Processing
            result.chief_history_num = result.combined_history(:, 1:6);
            result.deputy_history_num = result.combined_history(:, 7:12);
            result.absolute_state_history = util.ECI2RTN_history(result.deputy_history_num, result.chief_history_num);
            result.deputy_state_history_eci = result.deputy_history_num;
        end

        function statedot = wrapper_state_to_stateDot(self, t, state, result)
            chief_state_history = result.chief_history_num;
            t_vec = result.t_num;
            control_check = false;            

            if (t == 0) || (t >= t_vec(self.current_step + 1) && self.previous_time < t_vec(self.current_step + 1))
                if t == 0
                    [estimated_state, self.P] = our_algorithms.state_estimation(state, state, chief_state_history(1, :), self.dt, self.P, self.simulation_settings);
                else
                    [estimated_state, self.P] = our_algorithms.state_estimation(self.estimated_state_history(self.current_step-1, :)', state, chief_state_history(self.current_step-1, :), self.dt, self.P, self.simulation_settings);
                end
                self.estimated_state_history(self.current_step, :) = estimated_state;
                self.covariance_history(self.current_step, :, :) = self.P;


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CONTROL ALG
                chief_eci = chief_state_history(self.current_step+1, :)';
                deputy_rtn = util.ECI2RTN(state, chief_eci);
                result.dv = [result.dv; self.delta_v_tracker*ones(1,3)];
                burn_dv = 0.000001;
                old_v = deputy_rtn(4:6);

                deputy_rtn_new = deputy_rtn;
                deputy_rtn_new(4:6) = burn_dv * -deputy_rtn(1:3) / norm(deputy_rtn(1:3));
                self.delta_v_tracker = self.delta_v_tracker + norm(deputy_rtn(4:6) - old_v);

                old_eci = state;
                new_eci = util.RTN2ECI(deputy_rtn_new, chief_eci);
                eci_diff = new_eci - old_eci;
                v_diff = eci_diff(4:6);
                self.acceleration = v_diff / self.dt;

                self.last_control_time = t;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if t ~= t_vec(end)
                    self.current_step = self.current_step + 1;
                end
                self.previous_time = t;             
            end
            
            statedot = dynamics.two_body_dynamics_control(t, state, self.simulation_settings, self.acceleration);
        end
    end
end