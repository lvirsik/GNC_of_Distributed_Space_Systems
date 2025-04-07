function plotter(result, graphics_settings)
    if graphics_settings.orbit_eci
        plot_orbit_eci(result)
    end

    if graphics_settings.compare_numerical_vs_kepler
        compare_numerical_vs_kepler(result)
    end

    if graphics_settings.plot_orbital_elements.base_elems
        plot_orbital_elements(result, graphics_settings)
    end

    if graphics_settings.plot_eccentricity_vector
        plot_eccentricity_vector(result)
    end

    if graphics_settings.plot_ang_momentum_vector
        plot_ang_momentum_vector(result)
    end

    if graphics_settings.plot_specific_energy
        plot_specific_energy(result)
    end
end

function compare_numerical_vs_kepler(result)
    t = result.t_num;
    state_history_num = result.state_history_num;
    t_span = result.t_kep;
    state_history_kep = result.state_history_kep;

    state_history_kep = interp1(t_span, state_history_kep, t, 'linear');

    error_rtn = zeros(size(state_history_kep));
    for j = 1:size(state_history_kep, 1)
        error_rtn(j, :) = util.ECI2RTN(state_history_kep(j, :)) - util.ECI2RTN(state_history_num(j, :));
    end
    
    figure;
    hold on;
    plot(t, error_rtn)
    xlabel('Time (s)');
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

    figure;

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
        J2_analytical = zeros(length(result.state_history_num), 2);
        for j = 1:length(result.state_history_num)
            K = -(3 * constants.J2 * sqrt(constants.mu) * constants.earth_radius^2) / (2 * ((1 - e(j)^2)^2) * (a(j) ^ (7/2)));
            RAAN_j2 = result.initial_conditions(4) - (K * t(j) * cos(i(j)));
            w_j2 = result.initial_conditions(5) + ((K / 2) * t(j) * ((5 * (cos(i(j)) ^ 2)) - 1));
            J2_analytical(j, :) = [RAAN_j2, w_j2];
        end
    end

    figure;
    tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t, a, 'r'); grid on;
    xlabel('Time'); ylabel('a (km)');
    title('Semi-Major Axis');

    nexttile;
    plot(t, e, 'g'); grid on;
    xlabel('Time'); ylabel('e');
    title('Eccentricity');

    nexttile;
    plot(t, rad2deg(i), 'y'); grid on;
    xlabel('Time'); ylabel('Inclination (deg)');
    title('Inclination');

    nexttile;
    hold on;
    plot(t, rad2deg(RAAN), 'm'); grid on;
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t, rad2deg(J2_analytical(:,1)), 'k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    xlabel('Time'); ylabel('RAAN (deg)');
    title('RAAN');

    nexttile;
    hold on;
    plot(t, rad2deg(omega), 'c'); grid on;
    if graphics_settings.plot_orbital_elements.J2_mean_analytical_solution
        plot(t, rad2deg(J2_analytical(:,2)), 'k');
        legend('Numerical', 'J2 Mean Analytical');
    end
    xlabel('Time'); ylabel('\omega (deg)');
    title('Argument of Periapsis');
    

    nexttile;
    plot(t, rad2deg(nu), 'k'); grid on;
    xlabel('Time'); ylabel('\nu (deg)');
    title('True Anomaly');
end

function plot_eccentricity_vector(result)
    figure;
    hold on
    plot(result.t_num, result.ecc_vector_history(:,1), 'r'); 
    plot(result.t_num, result.ecc_vector_history(:,2), 'g');
    plot(result.t_num, result.ecc_vector_history(:,3), 'b');
    grid on;
    xlabel('Time (s)'); ylabel('Eccentricity Vector');
    title('Eccentricity Vector over Time');
    legend('X_comp', 'Y_comp', 'Z_comp');
end

function plot_ang_momentum_vector(result)
    figure;
    hold on
    plot(result.t_num, result.ang_mom_history(:,1), 'r'); 
    plot(result.t_num, result.ang_mom_history(:,2), 'g');
    plot(result.t_num, result.ang_mom_history(:,3), 'b');
    grid on;
    xlabel('Time (s)'); ylabel('Angular Momentum Vector');
    title('Angular Momentum Vector over Time');
    legend('X_comp', 'Y_comp', 'Z_comp');
end

function plot_specific_energy(result)
    figure;
    hold on
    plot(result.t_num, result.energy_history, 'b'); 
    grid on;
    xlabel('Time (s)'); ylabel('Specific Energy (J)');
    title('Specific Energy over Time');
end