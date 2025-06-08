classdef our_algorithms
    properties
        estimated_state_history = zeros(1, 6)
    end
    methods ( Static )
        function [estimated_state, covariance] = state_estimation(estimated_state, truth_state, chief_state, dt, P, simulation_settings, control_input)
        
            % Measurements: ECI Position, ECI Velocity, RTN Position, RTN Velocity
            measurement = our_algorithms.sensor_measurements(truth_state, chief_state);
            
            % Kalman Filter
            [estimated_state, covariance] = our_algorithms.kalman_filter(measurement, estimated_state, chief_state, dt, P, simulation_settings, control_input);
        end

        function [x_1, P_1] = kalman_filter(measurement, previous_state, chief_state, dt, P, simulation_settings, control_input)

            Q = diag([1e1, 1e1, 1e1, 1e3, 1e3, 1e3]);
            R = eye(6);

            % Measurements: ECI Position, ECI Velocity, RTN Position, RTN Velocity
            H = [[1,0,0,0,0,0];
                [0,1,0,0,0,0];
                [0,0,1,0,0,0];
                [0,0,0,1,0,0];
                [0,0,0,0,1,0];
                [0,0,0,0,0,1];
            ];

            % Predict Step
            A = our_algorithms.compute_linearized_dynamics(previous_state, simulation_settings);
            B = our_algorithms.compute_linearized_control(previous_state);
            F = eye(6) +  A * dt + B * u;
            x_0 = F * previous_state;
            P_0 = F * P * F' + Q;

            % Update Step
            z = H * x_0;
            y = measurement - z;
            K = P_0 * H' * inv((H * P_0 * H') + R);
            x_1 = x_0 + K * y;
            P_1 = (eye(6) - (K * H)) * P_0;
            
        end

        function measurement = sensor_measurements(truth_state, chief_state)
            pos_noise_std = 1.0;  % meters
            vel_noise_std = 0.1;  % meters/second
            pos_noise = pos_noise_std * randn(1, 3)';
            vel_noise = vel_noise_std * randn(1, 3)';

            measurement = truth_state;
            measurement(1:3) = truth_state(1:3) + pos_noise;  % add noise to position
            measurement(4:6) = truth_state(4:6) + vel_noise;  % add noise to velocity
        end

        function A = compute_linearized_dynamics(truth_state, simulation_settings)
            A = zeros(6, 6);
            for i=1:6
                state_plus = truth_state;
                state_minus = truth_state;

                h = 0.1;
                state_plus(i) = state_plus(i) + h;
                state_minus(i) = state_minus(i) - h;

                statedot_plus = dynamics.two_body_dynamics(0, state_plus, simulation_settings);
                statedot_minus = dynamics.two_body_dynamics(0, state_minus, simulation_settings);

                A(:,i) = ((statedot_plus - statedot_minus) / (2 * h));
            end
        end

        function B = compute_linearized_control(truth_state)
            B = zeros(6, 3);
            B(4, 1) = 1;
            B(5, 2) = 1;
            B(6, 3) = 1;
        end
    end
end