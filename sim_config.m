% Input initial conditions in orbital elments. Units of meters and degrees
initial_conditions_chief = [6780000, 0.0006, 51.6, 0, 0, 0]; % OEs
initial_conditions_chief = [6780000, 0, 0, 0, 0, 0]; 
initial_conditions_deputy = [0; 0; 0; 1; 1; 0]; % Position and velocity relative to chief in RTN

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
graphics_settings.plot_relative_position_deputy_comparison = true;
graphics_settings.plot_relative_position_error = true;

graphics_settings.plot_orbital_elements = struct();
graphics_settings.plot_orbital_elements.base_elems = false;
graphics_settings.plot_orbital_elements.J2_mean_analytical_solution = false;

graphics_settings.plot_eccentricity_vector = false;
graphics_settings.plot_ang_momentum_vector = false;
graphics_settings.plot_specific_energy = false;

num_orbits = 3;
time_step = 1;
disp('Running Case 1: Equal semi-major axis');
sim1 = simulator(initial_conditions_chief, initial_conditions_deputy, num_orbits, time_step, simulation_settings, graphics_settings);
sim1.run_simulator();

% Case 2: Create deputy with non-zero difference in semi-major axis
% This is done by adding a radial velocity component, which changes the energy
disp('Running Case 2: Non-zero difference in semi-major axis');
initial_conditions_deputy_case2 = [0; 0; 0; 10; 1; 0];
sim2 = simulator(initial_conditions_chief, initial_conditions_deputy_case2, num_orbits, time_step, simulation_settings, graphics_settings);
sim2.run_simulator();