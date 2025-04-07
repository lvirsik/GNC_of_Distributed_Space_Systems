function statedot = dynamics(t, state, simulation_settings)
% Take in state vector and return statedot vector (accelerations)
r = state(1:3);
v = state(4:6);

a = -constants.mu * r / (norm(r) ^ 3);

if simulation_settings.J2
    a_j2_x = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (1 - (5 * (r(3) / norm(r)) ^ 2)) * r(1);
    a_j2_y = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (1 - (5 * (r(3) / norm(r)) ^ 2)) * r(2);
    a_j2_z = 1.5 * constants.J2 * constants.mu * (constants.earth_radius ^ 2 / norm(r) ^ 5) * (3 - (5 * (r(3) / norm(r)) ^ 2)) * r(3);
    a_j2 = [a_j2_x; a_j2_y; a_j2_z];
    a = a + a_j2;
end

r_dot = v;
v_dot = a;

statedot = [r_dot; v_dot];

end