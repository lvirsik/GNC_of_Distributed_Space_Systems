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
                if self.simulation_settings.deputy
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
            plotter(result, self.graphics_settings)
        end
    end
end