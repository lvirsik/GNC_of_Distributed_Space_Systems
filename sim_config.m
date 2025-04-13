% Input initial conditions in orbital elments. Units of meters and degrees
initial_conditions_chief = [6780000, 0.0006, 51.6, 0, 0, 0]; % OEs
initial_conditions_chief = [6780000, 0, 0, 0, 0, 0]; 
initial_conditions_deputy = [0; 0; 0; 1; 1; 0]; % Position and velocity relative to chief in RTN
num_orbits = 3;
time_step = 1;

% Simulation settings
simulation_settings.numerical_propogation = true;
simulation_settings.keplerian_propogation = false;
simulation_settings.J2 = false;

simulation_settings.relative_deputy = true;
simulation_settings.absolute_deputy = true;


% Set graphics settings for graph output
graphics_settings.orbit_eci = false;
graphics_settings.compare_numerical_vs_kepler = false;

graphics_settings.plot_relative_position_deputy = true;

graphics_settings.plot_orbital_elements = struct();
graphics_settings.plot_orbital_elements.base_elems = false;
graphics_settings.plot_orbital_elements.J2_mean_analytical_solution = false;

graphics_settings.plot_eccentricity_vector = false;
graphics_settings.plot_ang_momentum_vector = false;
graphics_settings.plot_specific_energy = false;


sim = simulator(initial_conditions_chief, initial_conditions_deputy, num_orbits, time_step, simulation_settings, graphics_settings);
sim.run_simulator();