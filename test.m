% Test and visualize ECI <-> RTN round-trip conversion
clc; clear;

% Generate random test states
num_tests = 100;
eci_errors = zeros(num_tests, 1);

for i = 1:num_tests
    % Simulate random deputy and chief states (scaled for realism)
    r_chief = randn(3,1) * 7000e3;   % ~LEO radius
    v_chief = randn(3,1) * 7.5e3;    % ~LEO velocity
    chief_state = [r_chief; v_chief];

    r_deputy = r_chief + randn(3,1) * 1000;   % Offset by ~1 km
    v_deputy = v_chief + randn(3,1) * 1;      % Slight velocity offset
    deputy_state = [r_deputy; v_deputy];

    % Convert ECI -> RTN -> ECI
    state_rtn = util.ECI2RTN(deputy_state, chief_state);
    recovered_eci = util.RTN2ECI(state_rtn, chief_state);

    % Compute Euclidean norm of error
    eci_errors(i) = norm(recovered_eci - deputy_state);
end

% Plot the error
figure;
plot(1:num_tests, eci_errors * 1e3, 'b-o', 'LineWidth', 1.5);
xlabel('Test Index');
ylabel('Round-trip Error (mm)');
title('ECI → RTN → ECI Round-trip Error');
grid on;

% Display max and average error
fprintf('Max Error: %.6f mm\n', max(eci_errors) * 1e3);
fprintf('Mean Error: %.6f mm\n', mean(eci_errors) * 1e3);
