function plotter(result, graphics_settings)
    % Create figures directory in the current working directory
    figures_dir = 'figures';
    
    % Get the case identifier
    if isfield(result, 'maneuver_delta_v')
        case_str = 'Case3_';
    elseif isfield(result, 'deputy_in_rtn') && max(abs(result.deputy_in_rtn(:,1))) > 5000
        case_str = 'Case2_';
    else
        case_str = 'Case1_';
    end
    
    % Counter for figures
    fig_count = 0;
    
    % Close any existing figures to avoid interference
    close all;
    
    % Display all figures and save them
    if graphics_settings.orbit_eci
        fig_count = fig_count + 1;
        fig = figure('Name', 'OrbitECI');
        plot_orbit_eci(result);
        saveas(fig, ['./figures/', case_str, 'OrbitECI.png']);
    end

    if graphics_settings.compare_numerical_vs_kepler
        fig_count = fig_count + 1;
        fig = figure('Name', 'NumVsKep');
        compare_numerical_vs_kepler(result);
        saveas(fig, ['./figures/', case_str, 'NumVsKep.png']);
    end

    if graphics_settings.plot_relative_position_deputy
        % Position figure
        fig_count = fig_count + 1;
        fig_traj = figure('Name', 'RelativeTrajectory');
        
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,1); T = rho(:,2); N = rho(:,3);
        
        hold on;
        plot3(T, N, R, 'b', 'LineWidth', 2);
        xlabel('Tangential (m)');
        ylabel('Normal (m)');
        zlabel('Radial (m)');
        title('Deputy Trajectory in RTN Frame');
        grid on;
        view(3);
        
        saveas(fig_traj, ['./figures/', case_str, 'RelativeTrajectory.png']);
        
        % Velocity figure
        fig_count = fig_count + 1;
        fig_vel = figure('Name', 'RelativeVelocity');
        
        Rv = rho(:,4); Tv = rho(:,5); Nv = rho(:,6);
        
        hold on;
        plot3(Tv*100, Nv*100, Rv*100, 'r', 'LineWidth', 2);
        xlabel('Tangential (cm/s)');
        ylabel('Normal (cm/s)');
        zlabel('Radial (cm/s)');
        title('Deputy Velocity in RTN Frame');
        grid on;
        view(3);
        
        saveas(fig_vel, ['./figures/', case_str, 'RelativeVelocity.png']);
    end
    
    if isfield(graphics_settings, 'plot_relative_position_deputy_comparison') && graphics_settings.plot_relative_position_deputy_comparison && isfield(result, 'deputy_in_rtn')
        fig_count = fig_count + 1;
        fig = figure('Name', 'MethodsComparison');
        
        % For (2b) - relative motion from nonlinear equations
        rho = result.relative_state_history(:, 1:3);
        
        % For (2c) - relative motion computed from differencing individual orbits
        rho_diff = result.deputy_in_rtn(:, 1:3);
        
        hold on;
        plot3(rho(:,2), rho(:,3), rho(:,1), 'b', 'LineWidth', 2);
        plot3(rho_diff(:,2), rho_diff(:,3), rho_diff(:,1), 'r--', 'LineWidth', 1.5);
        xlabel('Tangential (m)');
        ylabel('Normal (m)');
        zlabel('Radial (m)');
        title('Deputy Trajectory in RTN Frame - Methods Comparison');
        legend('Nonlinear Relative Equations', 'Differencing Individual Orbits');
        grid on;
        view(3);
        
        saveas(fig, ['./figures/', case_str, 'MethodsComparison.png']);
    end
    
    if isfield(graphics_settings, 'plot_relative_position_error') && graphics_settings.plot_relative_position_error && isfield(result, 'deputy_in_rtn')
        fig_count = fig_count + 1;
        fig_error = figure('Name', 'MethodsError');
        
        % Calculate error between nonlinear equations and differencing orbits
        rho_nonlinear = result.relative_state_history(:, 1:3);
        rho_differencing = result.deputy_in_rtn(:, 1:3);
        
        % Ensure they're the same length for comparison
        min_length = min(size(rho_nonlinear, 1), size(rho_differencing, 1));
        rho_nonlinear = rho_nonlinear(1:min_length, :);
        rho_differencing = rho_differencing(1:min_length, :);
        
        % Calculate error in each component
        error = rho_nonlinear - rho_differencing;
        
        % Calculate absolute error
        abs_error = sqrt(sum(error.^2, 2));
        
        % Plot errors - using the current figure (fig_error)
        subplot(2,1,1);
        hold on;
        plot(result.t_num(1:min_length) / (60*60), error(:,1), 'r');
        plot(result.t_num(1:min_length) / (60*60), error(:,2), 'g');
        plot(result.t_num(1:min_length) / (60*60), error(:,3), 'b');
        grid on;
        xlabel('Time (hrs)');
        ylabel('Error (m)');
        title('Component Error Between Nonlinear Equations and Differencing Orbits');
        legend('Radial Error', 'Tangential Error', 'Normal Error');
        
        subplot(2,1,2);
        plot(result.t_num(1:min_length) / (60*60), abs_error, 'k');
        grid on;
        xlabel('Time (hrs)');
        ylabel('Absolute Error (m)');
        title('Absolute Error Between Methods');
        
        % Display statistics
        max_error = max(abs_error);
        mean_error = mean(abs_error);
        disp(['Maximum absolute error: ', num2str(max_error), ' meters']);
        disp(['Mean absolute error: ', num2str(mean_error), ' meters']);
        
        saveas(fig_error, ['./figures/', case_str, 'MethodsError.png']);
    end
    
    if isfield(graphics_settings, 'plot_maneuver_comparison') && graphics_settings.plot_maneuver_comparison && isfield(result, 'deputy_in_rtn_combined')
        fig_count = fig_count + 1;
        fig = figure('Name', 'ManeuverComparison');
        
        % Check that the necessary fields exist
        if ~isfield(result, 'deputy_in_rtn_combined')
            disp('No maneuver applied or combined trajectory not available');
        else
            % Original trajectory (pre-maneuver)
            rho_original = result.deputy_in_rtn(:, 1:3);
            
            % Combined trajectory (with maneuver)
            rho_combined = result.deputy_in_rtn_combined(:, 1:3);
            
            % Find maneuver point
            maneuver_idx = result.maneuver_time_index;
            maneuver_point = rho_original(maneuver_idx, :);
            
            hold on;
            plot3(rho_original(:,2), rho_original(:,3), rho_original(:,1), 'r', 'LineWidth', 1.5);
            plot3(rho_combined(:,2), rho_combined(:,3), rho_combined(:,1), 'b', 'LineWidth', 1.5);
            plot3(maneuver_point(2), maneuver_point(3), maneuver_point(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            
            xlabel('Tangential (m)');
            ylabel('Normal (m)');
            zlabel('Radial (m)');
            title('Deputy Trajectory Before and After Maneuver');
            legend('Original Trajectory (Unbounded)', 'Post-Maneuver Trajectory (Bounded)', 'Maneuver Location');
            grid on;
            view(3);
            
            % Display maneuver information
            annotation('textbox', [0.15, 0.05, 0.7, 0.1], 'String', ...
                {['Maneuver at t = ', num2str(result.maneuver_time/3600, '%.2f'), ' hours'], ...
                 ['Delta-v magnitude = ', num2str(result.maneuver_delta_v, '%.4f'), ' m/s']}, ...
                'FitBoxToText', 'on', 'BackgroundColor', 'white');
        end
        
        saveas(fig, ['./figures/', case_str, 'ManeuverComparison.png']);
    end

    if graphics_settings.plot_orbital_elements.base_elems
        fig_count = fig_count + 1;
        fig = figure('Name', 'OrbitalElements');
        plot_orbital_elements(result, graphics_settings);
        saveas(fig, ['./figures/', case_str, 'OrbitalElements.png']);
    end

    if graphics_settings.plot_eccentricity_vector
        fig_count = fig_count + 1;
        fig = figure('Name', 'EccentricityVector');
        plot_eccentricity_vector(result);
        saveas(fig, ['./figures/', case_str, 'EccentricityVector.png']);
    end

    if graphics_settings.plot_ang_momentum_vector
        fig_count = fig_count + 1;
        fig = figure('Name', 'AngularMomentum');
        plot_ang_momentum_vector(result);
        saveas(fig, ['./figures/', case_str, 'AngularMomentum.png']);
    end

    if graphics_settings.plot_specific_energy
        fig_count = fig_count + 1;
        fig = figure('Name', 'SpecificEnergy');
        plot_specific_energy(result);
        saveas(fig, ['./figures/', case_str, 'SpecificEnergy.png']);
    end
    
    disp(['Generated ' num2str(fig_count) ' figures and saved them to ./figures/ directory']);
end

function compare_numerical_vs_kepler(result)
    t = result.t_num;
    state_history_num = result.state_history_num;
    t_span = result.t_kep;
    state_history_kep = result.state_history_kep;

    state_history_kep = interp1(t_span, state_history_kep, t, 'linear');

    error_rtn = zeros(size(state_history_kep));
    for j = 1:size(state_history_kep, 1)
        R_eci2rtn = util.R_ECI2RTN(state_history_kep(j, :));
        r_rtn_kep = R_eci2rtn * state_history_kep(j, 1:3)';
        r_rtn_num = R_eci2rtn * state_history_num(j, 1:3)';


        f_dot = norm(cross(state_history_kep(j, 1:3), state_history_kep(j, 4:6))) / norm(state_history_kep(j, 1:3))^2;
        w = [0,0,f_dot]';
        v_rtn_kep = (R_eci2rtn * state_history_kep(j, 4:6)') - cross(w, r_rtn_kep);
        v_rtn_num = (R_eci2rtn * state_history_num(j, 4:6)') - cross(w, r_rtn_kep);

        error_rtn(j, :) = [r_rtn_kep; v_rtn_kep] - [r_rtn_num; v_rtn_num];
    end
    
    hold on;
    plot(t/(60*60), error_rtn)
    xlabel('Time (hrs)');
    ylabel('Absolute Error (m and m/s)');
    title('Absolute Error vs Time, RTN Frame');
    legend('Error R_r', 'Error R_t', 'Error R_n', 'Error V_r', 'Error V_t', 'Error V_n');
end

function plot_orbit_eci(result)
    state_history = result.state_history_num;

    % Extract positions
    x = state_history(:, 1);
    y = state_history(:, 2);
    z = state_history(:, 3);

    % Show Earth
    [xe, ye, ze] = sphere(50);
    xe = constants.earth_radius * xe;
    ye = constants.earth_radius * ye;
    ze = constants.earth_radius * ze;

    % Plot Earth
    surf(xe, ye, ze, 'FaceColor', 'blue', 'EdgeColor', 'none');
    hold on;

    % Plot Orbit
    plot3(x, y, z, 'c-', 'LineWidth', 2); 
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('ECI Trajectory');
    view(3);
    axis equal;
end

function plot_orbital_elements(result, graphics_settings)
    t = result.t_num;
    oe_history = result.oe_history;

    a     = oe_history(:, 1);
    e     = oe_history(:, 2);
    i     = oe_history(:, 3);
    RAAN  = oe_history(:, 4);  
    omega = oe_history(:, 5);
    nu    = oe_history(:, 6);  

    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        J2_analytical = zeros(length(result.state_history_num), 6);
        for j = 1:length(result.state_history_num)
            K = -(3 * constants.J2 * sqrt(constants.mu) * constants.earth_radius^2) / (2 * ((1 - e(j)^2)^2) * (a(j) ^ (7/2)));
            a_j2 = result.initial_conditions(1);
            e_j2 = result.initial_conditions(2);
            i_j2 = result.initial_conditions(3);
            RAAN_j2 = result.initial_conditions(4) - (K * t(j) * cos(i(j)));
            w_j2 = result.initial_conditions(5) + ((K / 2) * t(j) * ((5 * (cos(i(j)) ^ 2)) - 1)); r = util.OE2ECI(a(j), e(j), i(j), RAAN(j), omega(j), nu(j));
            v_j2 = result.initial_conditions(5) + (result.dt * ((sqrt(constants.mu * a(j)*(1 - e(j)^2)) / norm(r(1:3))^2) - ((K/2) * (5*cos(i(j))^2 - 1))));
            J2_analytical(j, :) = [a_j2, e_j2, i_j2, RAAN_j2, w_j2, v_j2];
        end
    end

    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    hold on;
    plot(t/(60*60), a, 'r'); grid on;
    xlabel('Time (hrs)'); ylabel('a (km)');
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), J2_analytical(:,1), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    title('Semi-Major Axis');

    nexttile;
    hold on;
    plot(t/(60*60), e, 'g'); grid on;
    xlabel('Time (hrs)'); ylabel('e');
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), J2_analytical(:,2), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    title('Eccentricity');

    nexttile;
    hold on;
    plot(t/(60*60), rad2deg(i), 'y'); grid on;
    xlabel('Time (hrs)'); ylabel('Inclination (deg)');
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), J2_analytical(:,3), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    title('Inclination');

    nexttile;
    hold on;
    plot(t/(60*60), rad2deg(RAAN), 'm'); grid on;
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), rad2deg(J2_analytical(:,4)), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    xlabel('Time (hrs)'); ylabel('RAAN (deg)');
    title('RAAN');

    nexttile;
    hold on;
    plot(t/(60*60), rad2deg(omega), 'c'); grid on;
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), rad2deg(J2_analytical(:,5)), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    xlabel('Time (hrs)'); ylabel('\omega (deg)');
    title('Argument of Periapsis');
    

    nexttile;
    hold on;
    plot(t/(60*60), rad2deg(nu), 'k'); grid on;
    xlabel('Time (hrs)'); ylabel('\nu (deg)');
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t/(60*60), rad2deg(J2_analytical(:,6)), '--k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    title('True Anomaly');
end

function plot_eccentricity_vector(result)
    hold on
    plot(result.t_num / (60*60), result.ecc_vector_history(:,1), 'r'); 
    plot(result.t_num / (60*60), result.ecc_vector_history(:,2), 'g');
    plot(result.t_num / (60*60), result.ecc_vector_history(:,3), 'b');
    grid on;
    xlabel('Time (hrs)'); ylabel('Eccentricity Vector');
    title('Eccentricity Vector over Time');
    legend('X_comp', 'Y_comp', 'Z_comp');
end

function plot_ang_momentum_vector(result)
    hold on
    plot(result.t_num / (60*60), result.ang_mom_history(:,1), 'r'); 
    plot(result.t_num / (60*60), result.ang_mom_history(:,2), 'g');
    plot(result.t_num / (60*60), result.ang_mom_history(:,3), 'b');
    grid on;
    xlabel('Time (hrs)'); ylabel('Angular Momentum Vector');
    title('Angular Momentum Vector over Time');
    legend('X_comp', 'Y_comp', 'Z_comp');
end

function plot_specific_energy(result)
    hold on
    plot(result.t_num / (60*60), result.energy_history, 'b'); 
    grid on;
    xlabel('Time (hrs)'); ylabel('Specific Energy (J)');
    title('Specific Energy over Time');
end