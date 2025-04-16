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
                
                if self.simulation_settings.absolute_deputy
                    [deputy_initial_state_eci, result.deputy_in_rtn, result.deputy_state_history, result.t_deputy] = util.propagate_deputy_orbit(self.initial_conditions_deputy, initial_state_eci, self.time_span, result.state_history_num, t, self.simulation_settings);
                                    
                    if self.simulation_settings.calculate_maneuver
                        [delta_v, optimal_time_index, maneuver_point] = util.calculate_drift_correction(result.deputy_state_history, result.state_history_num, result.t_num);
                        
                        result.maneuver_time = result.t_num(optimal_time_index);
                        result.maneuver_delta_v = delta_v;
                        result.maneuver_time_index = optimal_time_index;
                        result.maneuver_point = maneuver_point;
                        
                        if self.simulation_settings.apply_maneuver
                            [result.t_combined, result.deputy_state_combined, result.deputy_in_rtn_combined] = util.apply_maneuver(result, optimal_time_index, self.time_span, self.simulation_settings);
                        end
                    end
                end

                [result.oe_history, result.ecc_vector_history, result.ang_mom_history, result.energy_history] = util.calculate_orbit_history(result.state_history_num);
            end

            if self.simulation_settings.keplerian_propogation
                [result.state_history_kep, result.t_kep] = util.propagate_keplerian_orbit(a, e, i, RAAN, w, v, self.time_span, self.dt);
            end
            
            plotter(result, self.graphics_settings);
        end
    end
end