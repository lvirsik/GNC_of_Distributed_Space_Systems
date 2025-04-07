classdef simulator
    properties
        initial_conditions
        time_span
        simulation_settings
        graphics_settings
    end

    methods
        function self = simulator(initial_conditions, time_span, simulation_settings, graphics_settings)
            self.initial_conditions = initial_conditions;
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

            if self.simulation_settings.numerical_propogation
                initial_state = util.OE2ECI(a, e, i, RAAN, w, v);

                % Perform Simulation
                options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
                [t, state_history_num] = ode45(@(t, state_history_num) dynamics(t, state_history_num, self.simulation_settings), [0 self.time_span(end)], initial_state, options);
                result.state_history_num = state_history_num;
                result.t_num = t;

                % Gather orbital element history and other interesting information
                oe_history = zeros(length(state_history_num), 6);
                ecc_vector_history = zeros(length(state_history_num), 3);
                ang_mom_history = zeros(length(state_history_num), 3);
                energy_history = zeros(length(state_history_num), 1);
                for j = 1:length(state_history_num)
                    oe_history(j, :) = util.ECI2OE(state_history_num(j, :));
                    ecc_vector_history(j, :) = util.get_ecc_vector(state_history_num(j, :));
                    ang_mom_history(j, :) = util.get_ang_momentum(state_history_num(j, :));
                    energy_history(j, :) = util.get_energy(state_history_num(j, :));
                end
                result.oe_history = oe_history;
                result.ecc_vector_history = ecc_vector_history;
                result.ang_mom_history = ang_mom_history;
                result.energy_history = energy_history;
            end

            if self.simulation_settings.keplerian_propogation
                n = sqrt(constants.mu / a ^ 3);
                dt = self.time_span(2) - self.time_span(1);
                state_history_kep = zeros(length(self.time_span), 6);
                for j=1:length(self.time_span)
                    state_kep = util.OE2ECI(a, e, i, RAAN, w, v);
                    state_history_kep(j,:) = state_kep';
                    E = 2 * atan2(tan(v/2) * sqrt((1 - e) / (1 + e)), 1);
                    M = E - e*sin(E);
                    M_new = M + (n * dt);
                    E_new = util.MtoE(M_new, e, 10^(-6));
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