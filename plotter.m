function plotter(result, graphics_settings)
    % Display all figures and save them
    if graphics_settings.orbit_eci
        plot_orbit_eci(result);
    end

    if graphics_settings.compare_numerical_vs_kepler
        compare_numerical_vs_kepler(result);
    end

    if graphics_settings.plot_orbital_elements.base_elems
        plot_orbital_elements(result, graphics_settings);
    end

    if graphics_settings.plot_eccentricity_vector
        plot_eccentricity_vector(result);
    end

    if graphics_settings.plot_ang_momentum_vector
        plot_ang_momentum_vector(result);
    end

    if graphics_settings.plot_specific_energy
        plot_specific_energy(result);
    end

    if graphics_settings.plot_deputy.relative || graphics_settings.plot_deputy.absolute || graphics_settings.plot_deputy.hcw || graphics_settings.plot_deputy.ya || graphics_settings.plot_deputy.roe_eccentric
        plot_deputy(result, graphics_settings);
    end

    if graphics_settings.plot_deputy.manuvered
        plot_deputy_manuvered(result);
    end

    if graphics_settings.plot_quasi_oes
        plot_quasi_oes(result)
    end

    if graphics_settings.plot_quasi_roes
        plot_quasi_roes(result)
    end

    if graphics_settings.plot_quasi_both
        plot_quasi_both(result)
    end

    if graphics_settings.manuver_continuous
        plot_dragon_sim(result)
    end

    if graphics_settings.kalman_info
        kalman_info(result)
    end
end

function compare_numerical_vs_kepler(result)
    t = result.t_num;
    state_history_num = result.chief_history_num;
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
    state_history = result.chief_history_num;

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
        J2_analytical = zeros(length(result.chief_history_num), 6);
        for j = 1:length(result.chief_history_num)
            K = -(3 * constants.J2 * sqrt(constants.mu) * constants.earth_radius^2) / (2 * ((1 - e(j)^2)^2) * (a(j) ^ (7/2)));
            a_j2 = result.initial_conditions(1);
            e_j2 = result.initial_conditions(2);
            i_j2 = result.initial_conditions(3);
            RAAN_j2 = result.initial_conditions(4) - (K * t(j) * cos(i(j)));
            w_j2 = result.initial_conditions(5) + ((K / 2) * t(j) * ((5 * (cos(i(j)) ^ 2)) - 1)); r = util.OE2ECI([a(j), e(j), i(j), RAAN(j), omega(j), nu(j)]);
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

function plot_deputy(result, graphics_settings)
    % Position figure
    fig_traj = figure('Name', 'RelativeTrajectory');
    subplot(2,2,1);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,1); T = rho(:,2); N = rho(:,3);
        plot(T, R, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,1); T1 = rho(:,2); N1 = rho(:,3);
        plot(T1, R1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,1); T2 = rho(:,2); N2 = rho(:,3);
        plot(T2, R2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,1); T3 = rho(:,2); N3 = rho(:,3);
        plot(T3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,1); T4 = rho(:,2); N4 = rho(:,3);
        plot(T4, R4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(T5, R5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(T5, R5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Tangential (m)');
    ylabel('Radial (m)');
    legend(legend_names)
    grid on;

    subplot(2,2,2);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,1); T = rho(:,2); N = rho(:,3);
        plot(N, R, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,1); T1 = rho(:,2); N1 = rho(:,3);
        plot(N1, R1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,1); T2 = rho(:,2); N2 = rho(:,3);
        plot(N2, R2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,1); T3 = rho(:,2); N3 = rho(:,3);
        plot(N3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,1); T4 = rho(:,2); N4 = rho(:,3);
        plot(N4, R4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(N5, R5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(N5, R5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Normal (m)');
    ylabel('Radial (m)');
    legend(legend_names)
    grid on;

    subplot(2,2,3);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,1); T = rho(:,2); N = rho(:,3);
        plot(T, N, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,1); T1 = rho(:,2); N1 = rho(:,3);
        plot(T1, N1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,1); T2 = rho(:,2); N2 = rho(:,3);
        plot(T2, N2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,1); T3 = rho(:,2); N3 = rho(:,3);
        plot(T3, N3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,1); T4 = rho(:,2); N4 = rho(:,3);
        plot(T4, N4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(T5, N5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot(T5, N5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Tangential (m)');
    ylabel('Normal (m)');
    legend(legend_names)
    grid on;

    subplot(2,2,4);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,1); T = rho(:,2); N = rho(:,3);
        plot3(T, N, R, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,1); T1 = rho(:,2); N1 = rho(:,3);
        plot3(T1, N1, R1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,1); T2 = rho(:,2); N2 = rho(:,3);
        plot3(T2, N2, R2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,1); T3 = rho(:,2); N3 = rho(:,3);
        plot3(T3, N3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,1); T4 = rho(:,2); N4 = rho(:,3);
        plot3(T4, N4, R4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot3(T5, N5, R5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,1); T5 = rho(:,2); N5 = rho(:,3);
        plot3(T5, N5, R5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end

    xlabel('Tangential (m)');
    ylabel('Normal (m)');
    zlabel('Radial (m)');
    legend(legend_names)
    grid on;
    view(3);

    % Velocity figure
    fig_vel = figure('Name', 'RelativeVelocity');
    subplot(2,2,1);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,4); T = rho(:,5); N = rho(:,6);
        plot(T, R, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,4); T1 = rho(:,5); N1 = rho(:,6);
        plot(T1, R1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,4); T2 = rho(:,5); N2 = rho(:,6);
        plot(T2, R2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,4); T3 = rho(:,5); N3 = rho(:,6);
        plot(T3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,4); T4 = rho(:,5); N4 = rho(:,6);
        plot(T4, R4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(T5, R5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(T5, R5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Tangential Velocity (m/s)');
    ylabel('Radial Velocity (m/s)');
    legend(legend_names)
    grid on;

    subplot(2,2,2);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,4); T = rho(:,5); N = rho(:,6);
        plot(N, R, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,4); T1 = rho(:,5); N1 = rho(:,6);
        plot(N1, R1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,4); T2 = rho(:,5); N2 = rho(:,6);
        plot(N2, R2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,4); T3 = rho(:,5); N3 = rho(:,6);
        plot(N3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,4); T4 = rho(:,5); N4 = rho(:,6);
        plot(N4, R4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(N5, R5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(N5, R5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Normal Velocity (m/s)');
    ylabel('Radial Velocity (m/s)');
    legend(legend_names)
    grid on;

    subplot(2,2,3);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        R = rho(:,4); T = rho(:,5); N = rho(:,6);
        plot(T, N, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        R1 = rho(:,4); T1 = rho(:,5); N1 = rho(:,6);
        plot(T1, N1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        R2 = rho(:,4); T2 = rho(:,5); N2 = rho(:,6);
        plot(T2, N2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,4); T3 = rho(:,5); N3 = rho(:,6);
        plot(T3, N3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        R4 = rho(:,4); T4 = rho(:,5); N4 = rho(:,6);
        plot(T4, N4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(T5, N5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        R5 = rho(:,4); T5 = rho(:,5); N5 = rho(:,6);
        plot(T5, N5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Tangential Velocity (m/s)');
    ylabel('Normal Velocity (m/s)');
    legend(legend_names)
    grid on;

    subplot(2,2,4);
    hold on;
    legend_names = [];
    if graphics_settings.plot_deputy.relative
        rho = result.relative_state_history(:, 1:6);
        Rv = rho(:,4); Tv = rho(:,5); Nv = rho(:,6);
        plot3(Tv, Nv, Rv, 'b', 'LineWidth', 2);
        legend_names = [legend_names,"Relative Dynamics"];
    end
    if graphics_settings.plot_deputy.absolute
        rho = result.absolute_state_history(:, 1:6);
        Rv1 = rho(:,4); Tv1 = rho(:,5); Nv1 = rho(:,6);
        plot3(Tv1, Nv1, Rv1, 'r', 'LineWidth', 2);
        legend_names = [legend_names,"2Body Dynamics"];
    end
    if graphics_settings.plot_deputy.hcw
        rho = result.hcw_state_history(:, 1:6);
        Rv2 = rho(:,4); Tv2 = rho(:,5); Nv2 = rho(:,6);
        plot3(Tv2, Nv2, Rv2, 'g', 'LineWidth', 2);
        legend_names = [legend_names,"HCW"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        Rv3 = rho(:,4); Tv3 = rho(:,5); Nv3 = rho(:,6);
        plot3(Tv3, Nv3, Rv3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA"];
    end
    if graphics_settings.plot_deputy.roe_circular
        rho = result.roe_circular_state_history(:, 1:6);
        Rv4 = rho(:,4); Tv4 = rho(:,5); Nv4 = rho(:,6);
        plot3(Tv4, Nv4, Rv4, 'm', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Circular"];
    end
    if graphics_settings.plot_deputy.roe_eccentric
        rho = result.roe_eccentric_state_history(:, 1:6);
        Rv5 = rho(:,4); Tv5 = rho(:,5); Nv5 = rho(:,6);
        plot3(Tv5, Nv5, Rv5, 'y', 'LineWidth', 2);
        legend_names = [legend_names, "ROE Eccentric"];
    end
    if graphics_settings.plot_deputy.estimated_state
        rho = result.estimated_state_history(:, 1:6);
        Rv5 = rho(:,4); Tv5 = rho(:,5); Nv5 = rho(:,6);
        plot3(Tv5, Nv5, Rv5, 'b', 'LineWidth', 2);
        legend_names = [legend_names, "Estimated State"];
    end
    xlabel('Tangential (m/s)');
    ylabel('Normal (m/s)');
    zlabel('Radial (m/s)');
    title('Deputy Velocity in RTN Frame');
    grid on;
    legend(legend_names)
    view(3);
end

function plot_deputy_manuvered(result)
    % Position figure
    fig_traj = figure('Name', 'RelativeTrajectory');
    hold on;
    rho = result.deputy_manuvered(:, 1:6);
    R = rho(:,1); T = rho(:,2); N = rho(:,3);
    plot3(T, N, R, 'b', 'LineWidth', 2);
    xlabel('Tangential (m)');
    ylabel('Normal (m)');
    zlabel('Radial (m)');
    title('Deputy Trajectory in RTN Frame');
    grid on;
    view(3);

    % Velocity figure
    fig_vel = figure('Name', 'RelativeVelocity');
    hold on;
    rho = result.deputy_manuvered(:, 1:6);
    Rv = rho(:,4); Tv = rho(:,5); Nv = rho(:,6);
    plot3(Tv, Nv, Rv, 'b', 'LineWidth', 2);
    xlabel('Tangential (m/s)');
    ylabel('Normal (m/s)');
    zlabel('Radial (m/s)');
    title('Deputy Velocity in RTN Frame');
    grid on;
    view(3);
end

function plot_quasi_oes(result)
    % Transform deputy orbit to RTN frame relative to chief
    t = result.t_num;
    qOEs = zeros(length(result.t_num), 6);
    for i = 1:length(result.t_num)
        qOEs(i,:) = util.ECI2qOE(result.deputy_state_history_eci(i,:));
    end

    % Low-pass filter parameters
    sample_rate = 1 / mean(diff(t));       
    cutoff_freq = 1 / 11000; 
    [b, a] = butter(3, cutoff_freq / (0.5 * sample_rate), 'low');
    qOEs_mean = zeros(size(qOEs));
    for k = 1:2 
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end
    for k = 3:6
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end


    a     = qOEs(:, 1);     a_m = qOEs_mean(:, 1);
    l     = qOEs(:, 2);     l_m = qOEs_mean(:, 2);
    ex    = qOEs(:, 3);     ex_m = qOEs_mean(:, 3);
    ey    = qOEs(:, 4);     ey_m = qOEs_mean(:, 4);
    ix    = qOEs(:, 5);     ix_m = qOEs_mean(:, 5);
    iy    = qOEs(:, 6);     iy_m = qOEs_mean(:, 6);
    figure();
    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    t_hr = t / 3600;
    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile; hold on;
    plot(t_hr, a, 'r'); plot(t_hr, a_m, '--r');
    xlabel('Time (hrs)'); ylabel('a (km)'); title('Semi-Major Axis'); grid on;

    nexttile; hold on;
    plot(t_hr, l, 'g'); plot(t_hr, l_m, '--g');
    xlabel('Time (hrs)'); ylabel('l'); title('Mean Argument of Latitude'); grid on;

    nexttile; hold on;
    plot(t_hr, ex, 'y'); plot(t_hr, ex_m, '--y');
    xlabel('Time (hrs)'); ylabel('e_x'); title('Eccentricity X'); grid on;

    nexttile; hold on;
    plot(t_hr, ey, 'm'); plot(t_hr, ey_m, '--m');
    xlabel('Time (hrs)'); ylabel('e_y'); title('Eccentricity Y'); grid on;

    nexttile; hold on;
    plot(t_hr, ix, 'c'); plot(t_hr, ix_m, '--c');
    xlabel('Time (hrs)'); ylabel('i_x'); title('Inclination X'); grid on;

    nexttile; hold on;
    plot(t_hr, iy, 'k'); plot(t_hr, iy_m, '--k');
    xlabel('Time (hrs)'); ylabel('i_y'); title('Inclination Y'); grid on;
end

function plot_quasi_roes(result)
    % Transform deputy orbit to RTN frame relative to chief
    t = result.t_num;
    qOEs = zeros(length(result.t_num), 6);
    for i = 1:length(result.t_num)
        qOEs(i,:) = util.ECI2ROE(result.chief_history_num(i,:), result.deputy_state_history_eci(i,:));
    end

    % Low-pass filter parameters
    sample_rate = 1 / mean(diff(t));       
    cutoff_freq = 1 / 11000; 
    [b, a] = butter(3, cutoff_freq / (0.5 * sample_rate), 'low');
    qOEs_mean = zeros(size(qOEs));
    for k = 1:2 
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end
    for k = 3:6
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end


    a     = qOEs(:, 1);     a_m = qOEs_mean(:, 1);
    l     = qOEs(:, 2);     l_m = qOEs_mean(:, 2);
    ex    = qOEs(:, 3);     ex_m = qOEs_mean(:, 3);
    ey    = qOEs(:, 4);     ey_m = qOEs_mean(:, 4);
    ix    = qOEs(:, 5);     ix_m = qOEs_mean(:, 5);
    iy    = qOEs(:, 6);     iy_m = qOEs_mean(:, 6);
    figure();
    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    t_hr = t / 3600;
    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile; hold on;
    plot(t_hr, a, 'r'); plot(t_hr, a_m, '--r');
    xlabel('Time (hrs)'); ylabel('da'); title('da'); grid on;

    nexttile; hold on;
    plot(t_hr, l, 'g'); plot(t_hr, l_m, '--g');
    xlabel('Time (hrs)'); ylabel('dl'); title('dl'); grid on;

    nexttile; hold on;
    plot(t_hr, ex, 'y'); plot(t_hr, ex_m, '--y');
    xlabel('Time (hrs)'); ylabel('de_x'); title('dex'); grid on;

    nexttile; hold on;
    plot(t_hr, ey, 'm'); plot(t_hr, ey_m, '--m');
    xlabel('Time (hrs)'); ylabel('de_y'); title('dey'); grid on;

    nexttile; hold on;
    plot(t_hr, ix, 'c'); plot(t_hr, ix_m, '--c');
    xlabel('Time (hrs)'); ylabel('di_x'); title('dix'); grid on;

    nexttile; hold on;
    plot(t_hr, iy, 'k'); plot(t_hr, iy_m, '--k');
    xlabel('Time (hrs)'); ylabel('di_y'); title('diy'); grid on;
end

function plot_quasi_both(result)
    % Transform deputy orbit to RTN frame relative to chief
    t = result.t_num;
    qOEs = zeros(length(result.t_num), 6);
    for i = 1:length(result.t_num)
        qOEs(i,:) = util.ECI2ROE(result.chief_history_num(i,:), result.deputy_state_history_eci(i,:));
    end

    % Low-pass filter parameters
    sample_rate = 1 / mean(diff(t));       
    cutoff_freq = 1 / 11000; 
    [b, a] = butter(3, cutoff_freq / (0.5 * sample_rate), 'low');
    qOEs_mean = zeros(size(qOEs));
    for k = 1:2 
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end
    for k = 3:6
        qOEs_mean(:,k) = filtfilt(b, a, qOEs(:,k));
    end

    da     = qOEs(:, 1);     da_m = qOEs_mean(:, 1);
    dl     = qOEs(:, 2);     dl_m = qOEs_mean(:, 2);
    dex    = qOEs(:, 3);     dex_m = qOEs_mean(:, 3);
    dey    = qOEs(:, 4);     dey_m = qOEs_mean(:, 4);
    dix    = qOEs(:, 5);     dix_m = qOEs_mean(:, 5);
    diy    = qOEs(:, 6);     diy_m = qOEs_mean(:, 6);

    figure;
    tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    hold on
    plot(dex, dey, 'b');
    plot(dex_m, dey_m, '--b');
    xlabel('\Delta e_x'); ylabel('\Delta e_y');
    title('\Delta e_x vs \Delta e_y'); grid on; axis equal;

    nexttile;
    hold on
    plot(dix, diy, 'r');
    plot(dix_m, diy_m, '--r');
    xlabel('\Delta i_x'); ylabel('\Delta i_y');
    title('\Delta i_x vs \Delta i_y'); grid on; axis equal;

    nexttile;
    hold on
    plot(dl, da, 'k');
    plot(dl_m, da_m, '--k');
    xlabel('\Delta l'); ylabel('\Delta a');
    title('\Delta l vs \Delta a'); grid on; axis equal;

end

function plot_dragon_sim(result)
    figure;
    plot(result.t_num, result.dv);
    xlabel('Time (s)');
    ylabel('Delta V');
    title('Delta V vs Time');
    grid on;
end

function kalman_info(result)
    % Extract core data
    t = result.t_num;
    x_true = result.absolute_state_history;
    x_est = result.estimated_state_history;
    P_hist = result.covariance_history;
    N = length(t);
    
    figure;
    for i = 1:6
        subplot(3, 2, i);
        plot(t, x_true(:, i), 'k', 'LineWidth', 1.2); hold on;
        plot(t, x_est(:, i), 'b--', 'LineWidth', 1.2);
        xlabel('Time [s]');
        ylabel(sprintf('State %d', i));
        legend('True', 'Estimated');
        grid on;
    end
    sgtitle('True vs Estimated States');

    error = x_est - x_true;
    
    pos_sigma = zeros(N, 3);
    vel_sigma = zeros(N, 3);
    
    for i = 1:N
        P = squeeze(P_hist(i, :, :));
        pos_sigma(i, :) = sqrt(diag(P(1:3, 1:3)))';
        vel_sigma(i, :) = sqrt(diag(P(4:6, 4:6)))';
    end
    

    figure;
    for j = 1:3
        subplot(3,1,j);
        plot(t, error(:, j), 'b'); hold on;
        plot(t, pos_sigma(:, j), 'k--');
        plot(t, -pos_sigma(:, j), 'k--');
        xlabel('Time [s]');
        ylabel(sprintf('Position Error State %d [m]', j));
        title(sprintf('Position Error Component %d with ±1σ Bounds', j));
        grid on;
    end
    sgtitle('Position Errors with ±1σ Bounds');

    figure;
    for j = 1:3
        subplot(3,1,j);
        plot(t, error(:, j+3), 'r'); hold on;
        plot(t, vel_sigma(:, j), 'k--');
        plot(t, -vel_sigma(:, j), 'k--');
        xlabel('Time [s]');
        ylabel(sprintf('Velocity Error State %d [m/s]', j));
        title(sprintf('Velocity Error Component %d with ±1σ Bounds', j));
        grid on;
    end
    sgtitle('Velocity Errors with ±1σ Bounds');

    steady_indices = round(0.9 * N):N;
    mean_error = mean(error(steady_indices, :));
    std_error = std(error(steady_indices, :));

    fprintf('\nSteady-State Error Stats (last 10%% of simulation):\n');
    for i = 1:6
        fprintf('State %d: Mean = %.4f, Std = %.4f\n', i, mean_error(i), std_error(i));
    end

    if isfield(result, 'pre_fit_residuals') && isfield(result, 'post_fit_residuals')
        pre_fit = result.pre_fit_residuals;
        post_fit = result.post_fit_residuals;

        figure;
        subplot(2,1,1);
        plot(t, pre_fit, 'b'); hold on;
        if isfield(result, 'injected_noise_std')
            noise_std = result.injected_noise_std;
            if isscalar(noise_std)
                plot(t, noise_std * ones(size(t)), 'k--');
                plot(t, -noise_std * ones(size(t)), 'k--');
            else
                plot(t, noise_std, 'k--');
                plot(t, -noise_std, 'k--');
            end
        end
        xlabel('Time [s]'); ylabel('Pre-fit Residual');
        title('Pre-fit Residuals vs Injected Noise');
        grid on;

        subplot(2,1,2);
        plot(t, post_fit, 'r'); hold on;
        if isfield(result, 'injected_noise_std')
            if isscalar(noise_std)
                plot(t, noise_std * ones(size(t)), 'k--');
                plot(t, -noise_std * ones(size(t)), 'k--');
            else
                plot(t, noise_std, 'k--');
                plot(t, -noise_std, 'k--');
            end
        end
        xlabel('Time [s]'); ylabel('Post-fit Residual');
        title('Post-fit Residuals vs Injected Noise');
        grid on;

        sgtitle('Residuals Analysis');
    end
end
