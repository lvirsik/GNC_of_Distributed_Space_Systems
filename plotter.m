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

    if graphics_settings.plot_deputy.relative || graphics_settings.plot_deputy.absolute || graphics_settings.plot_deputy.hcw || graphics_settings.plot_deputy.ya
        plot_deputy(result, graphics_settings);
    end

    if graphics_settings.plot_deputy.manuvered
        plot_deputy_manuvered(result);
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
        legend_names = [legend_names,"HCW Dynamics"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        R3 = rho(:,1); T3 = rho(:,2); N3 = rho(:,3);
        plot3(T3, N3, R3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA Dynamics"];
    end

    xlabel('Tangential (m)');
    ylabel('Normal (m)');
    zlabel('Radial (m)');
    legend(legend_names)
    title('Deputy Trajectory in RTN Frame');
    grid on;
    view(3);

    % Velocity figure
    fig_vel = figure('Name', 'RelativeVelocity');
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
        legend_names = [legend_names,"HCW Dynamics"];
    end
    if graphics_settings.plot_deputy.ya
        rho = result.ya_state_history(:, 1:6);
        Rv3 = rho(:,4); Tv3 = rho(:,5); Nv3 = rho(:,6);
        plot3(Tv3, Nv3, Rv3, 'c', 'LineWidth', 2);
        legend_names = [legend_names,"YA Dynamics"];
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