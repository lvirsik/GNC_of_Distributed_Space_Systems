% Input initial conditions in orbital elments. Units of meters and degrees
initial_conditions_chief = [6780000, 0.0006, 51.6, 0, 0, 0]; % OEs
%initial_conditions_chief = [6780000, 0.001, 0.1, 0, 0, 0];
initial_conditions_deputy = [0; 0; 0; 5; 5; 1];
% initial_conditions_deputy = util.ROE2RTN([0; 100; 50; 100; 30; 200], initial_conditions_chief);

% Time step and number of orbits
num_orbits = 1;
time_step = 1;

% Simulation settings

simulation_settings.manuver_continuous = false;
simulation_settings.manuver_instant = true;

simulation_settings.numerical_propogation = true;
simulation_settings.keplerian_propogation = false;
simulation_settings.J2 = true;
simulation_settings.relative_deputy = false;
simulation_settings.absolute_deputy = true;
simulation_settings.hcw_deputy = false;
simulation_settings.ya_deputy = false;
simulation_settings.roe_circular_deputy  = false;
simulation_settings.roe_eccentric_deputy = false;

simulation_settings.create_bounded_motion = false;

% Set graphics settings
graphics_settings.manuver_continuous = false;
simulation_settings.manuver_instant = true;

graphics_settings.orbit_eci = false;
graphics_settings.compare_numerical_vs_kepler = false;

graphics_settings.plot_deputy = struct();
graphics_settings.plot_deputy.relative = false;
graphics_settings.plot_deputy.absolute = true;
graphics_settings.plot_deputy.hcw = false;
graphics_settings.plot_deputy.ya = false;
graphics_settings.plot_deputy.roe_circular  = false;
graphics_settings.plot_deputy.roe_eccentric = false;

graphics_settings.plot_quasi_oes = false;
graphics_settings.plot_quasi_roes = false;
graphics_settings.plot_quasi_both = false;

graphics_settings.plot_deputy.manuvered = false;

graphics_settings.plot_orbital_elements = struct();
graphics_settings.plot_orbital_elements.base_elems = false;
graphics_settings.plot_orbital_elements.J2_mean_analytical_solution = false;

graphics_settings.plot_eccentricity_vector = false;
graphics_settings.plot_ang_momentum_vector = false;
graphics_settings.plot_specific_energy = false;

sim = simulator(initial_conditions_chief, initial_conditions_deputy, num_orbits, time_step, simulation_settings, graphics_settings);
sim.run_simulator();