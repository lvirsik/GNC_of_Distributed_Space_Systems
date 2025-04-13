% Input initial conditions in orbital elments. Units of meters and degrees
initial_conditions = [6780000, 0.0006, 51.6, 0, 0, 0]; % OEs
num_orbits = 10;
time_step = 1;

% Simulation settings
simulation_settings.numerical_propogation = true;
simulation_settings.keplerian_propogation = false;
simulation_settings.J2 = true;

% Set graphics settings for graph output
graphics_settings.orbit_eci = true;
graphics_settings.compare_numerical_vs_kepler = false;

graphics_settings.plot_orbital_elements = struct();
graphics_settings.plot_orbital_elements.base_elems = false;
graphics_settings.plot_orbital_elements.J2_mean_analytical_solution = false;

graphics_settings.plot_eccentricity_vector = false;
graphics_settings.plot_ang_momentum_vector = false;
graphics_settings.plot_specific_energy = false;


sim = simulator(initial_conditions, num_orbits, time_step, simulation_settings, graphics_settings);
sim.run_simulator();