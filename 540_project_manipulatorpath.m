clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
gamma = 5;                % Adaptation gain for parameter updates
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target

% Thruster force magnitudes
forward_thrust_magnitude = 0.1; % Maximum forward thrust

% Gains
k_theta = 0.2;    % Gain for heading error correction
k_l = 0.5;        % Gain for lateral deviation correction
k_v = 1;          % Gain for velocity control

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_t, v, omega, q1, q2, dq1, dq2]
state = [0; 0; pi/4; 0.1; 0; deg2rad(30); deg2rad(-45); 0; 0];

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length);
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Manipulator Dynamics');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Define spacecraft and manipulator shapes
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]';
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);

% Main simulation loop
for i = 1:num_steps
    % Extract state variables
    s = state(1); l = state(2); theta_t = state(3);
    v = state(4); omega = state(5);
    q1 = state(6); q2 = state(7); dq1 = state(8); dq2 = state(9);

    % Scale s to index the path array
    index = 1 + (s / path_length) * (path_resolution - 1);
    index = max(1, min(round(index), path_resolution)); % Clamp to valid range

    % Retrieve current position and tangent
    current_position = path(:, index);
    tangent = path_tangent(:, index);

    % Compute the satellite's position in global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Control velocity and lateral deviation
    torque = -k_theta * theta_t; % Control heading error

    % Compute manipulator reaction torques
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];
    H = [m1_hat * l1^2 + m2_hat * (l1^2 + l2^2), m2_hat * l1 * l2 * cos(q2);
         m2_hat * l1 * l2 * cos(q2), m2_hat * l2^2];
    C = [0, -m2_hat * l1 * l2 * sin(q2) * dq2;
         m2_hat * l1 * l2 * sin(q2) * dq1, 0];
    tau_manipulator = -J' * (Kp_end_effector * [l; 0] + Kd_end_effector * [v; 0]); % End-effector torque
    
    % Debugging checks
    assert(all(size(H) == [2, 2]), 'H must be a 2x2 matrix!');
    assert(all(size(C) == [2, 2]), 'C must be a 2x2 matrix!');
    assert(all(size(tau_manipulator) == [2, 1]), 'tau_manipulator must be a 2x1 vector!');

    ddq = pinv(H) * (tau_manipulator - C * [dq1; dq2]); % Joint accelerations
    assert(all(size(ddq) == [2, 1]), 'ddq must be a 2x1 vector!');

    % Update joint states
    dq1 = dq1 + ddq(1) * dt;
    dq2 = dq2 + ddq(2) * dt;
    q1 = q1 + dq1 * dt;
    q2 = q2 + dq2 * dt;

    % Apply manipulator torque to satellite angular dynamics
    tau_reaction = -sum(tau_manipulator); % Reaction torque from manipulator
    angular_acceleration = (torque + tau_reaction) / inertia_satellite;

    % Update Frenet-Serret dynamics
    s_dot = v * cos(theta_t);
    l_dot = v * sin(theta_t);
    theta_t_dot = omega;
    v_dot = thrust / mass_satellite;
    omega_dot = angular_acceleration;

    % Update the state
    state = state + [s_dot; l_dot; theta_t_dot; v_dot; omega_dot; dq1; dq2; ddq] * dt;

    % Update visualization
    R = [cos(theta_t), -sin(theta_t); sin(theta_t), cos(theta_t)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
                     'YData', rotated_satellite(2, :) + position(2));
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    p2 = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];
    set(h_link1, 'XData', [position(1), p1(1)], 'YData', [position(2), p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);

    pause(0.05);
end
hold off;

%%

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
gamma = 5;                % Adaptation gain for parameter updates
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target
velocity_limit = 1.0;     % Maximum translational velocity
omega_limit = 1.0;        % Maximum angular velocity
damping_factor = 0.1;     % Damping factor for velocity and angular velocity

% Thruster force magnitudes
forward_thrust_magnitude = 0.1; % Maximum forward thrust

% Gains
k_theta = 0.2;    % Gain for heading error correction
k_l = 0.5;        % Gain for lateral deviation correction
k_v = 1;          % Gain for velocity control

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_t, v, omega, q1, q2, dq1, dq2]
state = [0; 0; pi/4; 0.1; 0; deg2rad(30); deg2rad(-45); 0; 0];

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length);
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Manipulator Dynamics');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Define spacecraft and manipulator shapes
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]';
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);

% Main simulation loop
for i = 1:num_steps
    % Extract state variables
    s = state(1); l = state(2); theta_t = state(3);
    v = state(4); omega = state(5);
    q1 = state(6); q2 = state(7); dq1 = state(8); dq2 = state(9);

    % Scale s to index the path array
    index = 1 + (s / path_length) * (path_resolution - 1);
    index = max(1, min(round(index), path_resolution)); % Clamp to valid range

    % Retrieve current position and tangent
    current_position = path(:, index);
    tangent = path_tangent(:, index);

    % Compute the satellite's position in global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Control velocity and lateral deviation
    torque = -k_theta * theta_t; % Control heading error

    % Compute manipulator reaction torques
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];
    H = [m1_hat * l1^2 + m2_hat * (l1^2 + l2^2), m2_hat * l1 * l2 * cos(q2);
         m2_hat * l1 * l2 * cos(q2), m2_hat * l2^2];
    C = [0, -m2_hat * l1 * l2 * sin(q2) * dq2;
         m2_hat * l1 * l2 * sin(q2) * dq1, 0];
    tau_manipulator = -J' * (Kp_end_effector * [l; 0] + Kd_end_effector * [v; 0]); % End-effector torque

    ddq = pinv(H) * (tau_manipulator - C * [dq1; dq2]); % Joint accelerations

    % Update joint states
    dq1 = dq1 + ddq(1) * dt;
    dq2 = dq2 + ddq(2) * dt;
    q1 = q1 + dq1 * dt;
    q2 = q2 + dq2 * dt;

    % Apply manipulator torque to satellite angular dynamics
    tau_reaction = -sum(tau_manipulator); % Reaction torque from manipulator
    angular_acceleration = (torque + tau_reaction) / inertia_satellite;

    % Update Frenet-Serret dynamics
    s_dot = v * cos(theta_t);
    l_dot = v * sin(theta_t);
    theta_t_dot = omega;
    v_dot = thrust / mass_satellite - damping_factor * v; % Add damping to velocity
    omega_dot = angular_acceleration - damping_factor * omega; % Add damping to angular velocity

    % Enforce velocity limits
    v = max(min(v + v_dot * dt, velocity_limit), -velocity_limit);
    omega = max(min(omega + omega_dot * dt, omega_limit), -omega_limit);

    % Update the state
    state = [s + s_dot * dt; l + l_dot * dt; theta_t + theta_t_dot * dt; ...
             v; omega; q1; q2; dq1; dq2];

    % Update visualization
    R = [cos(theta_t), -sin(theta_t); sin(theta_t), cos(theta_t)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
                     'YData', rotated_satellite(2, :) + position(2));
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    p2 = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];
    set(h_link1, 'XData', [position(1), p1(1)], 'YData', [position(2), p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);

    pause(0.05);
end
hold off;

%%

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
gamma = 5;                % Adaptation gain for parameter updates
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target
velocity_limit = 1.0;     % Maximum translational velocity
omega_limit = 1.0;        % Maximum angular velocity
dq_limit = 5;             % Maximum joint velocity
q_limit = pi;             % Joint angle limit
damping_factor = 0.1;     % Damping factor for velocity and angular velocity

% Thruster force magnitudes
forward_thrust_magnitude = 0.1; % Maximum forward thrust

% Gains
k_theta = 0.2;    % Gain for heading error correction
k_l = 0.5;        % Gain for lateral deviation correction
k_v = 1;          % Gain for velocity control

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_t, v, omega, q1, q2, dq1, dq2]
state = [0; 0; pi/4; 0.1; 0; deg2rad(30); deg2rad(-45); 0; 0];

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length);
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Manipulator Dynamics');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Define spacecraft and manipulator shapes
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]';
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5); % Trajectory plot

% Main simulation loop
for i = 1:num_steps
    % Extract state variables
    s = state(1); l = state(2); theta_t = state(3);
    v = state(4); omega = state(5);
    q1 = state(6); q2 = state(7); dq1 = state(8); dq2 = state(9);

    % Scale s to index the path array
    index = 1 + (s / path_length) * (path_resolution - 1);
    index = max(1, min(round(index), path_resolution)); % Clamp to valid range

    % Retrieve current position and tangent
    current_position = path(:, index);
    tangent = path_tangent(:, index);

    % Compute the satellite's position in global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;

    % Path error (lateral deviation)
    path_error_vector = position - current_position; % Vector from path point to spacecraft
    path_error = dot([-tangent(2); tangent(1)], path_error_vector); % Signed lateral deviation
    path_error_history(i) = path_error;

    % % Retrieve current position and tangent
    % current_position = path(:, index);
    % tangent = path_tangent(:, index);
    % 
    % % Compute the satellite's position in global frame
    % position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    % position_history(:, i) = position;

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Control velocity and lateral deviation
    torque = -k_theta * theta_t; % Control heading error

    % Compute manipulator reaction torques
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];
    H = [m1_hat * l1^2 + m2_hat * (l1^2 + l2^2), m2_hat * l1 * l2 * cos(q2);
         m2_hat * l1 * l2 * cos(q2), m2_hat * l2^2];
    C = [0, -m2_hat * l1 * l2 * sin(q2) * dq2;
         m2_hat * l1 * l2 * sin(q2) * dq1, 0];
    tau_manipulator = -J' * (Kp_end_effector * [l; 0] + Kd_end_effector * [v; 0]); % End-effector torque

    ddq = pinv(H) * (tau_manipulator - C * [dq1; dq2]); % Joint accelerations

    % Enforce limits on joint velocities and angles
    dq1 = max(min(dq1 + ddq(1) * dt, dq_limit), -dq_limit);
    dq2 = max(min(dq2 + ddq(2) * dt, dq_limit), -dq_limit);
    q1 = max(min(q1 + dq1 * dt, q_limit), -q_limit);
    q2 = max(min(q2 + dq2 * dt, q_limit), -q_limit);

    % Apply manipulator torque to satellite angular dynamics
    tau_reaction = sum(tau_manipulator); % Reaction torque from manipulator
    angular_acceleration = (torque + tau_reaction) / inertia_satellite;

    % Update Frenet-Serret dynamics
    s_dot = v * cos(theta_t);
    l_dot = v * sin(theta_t);
    theta_t_dot = omega;
    v_dot = thrust / mass_satellite - damping_factor * v; % Add damping to velocity
    omega_dot = angular_acceleration - damping_factor * omega; % Add damping to angular velocity

    % Enforce velocity limits
    v = max(min(v + v_dot * dt, velocity_limit), -velocity_limit);
    omega = max(min(omega + omega_dot * dt, omega_limit), -omega_limit);

    % Update the state
    state = [s + s_dot * dt; l + l_dot * dt; theta_t + theta_t_dot * dt; ...
             v; omega; q1; q2; dq1; dq2];

    % Update visualization
    R = [cos(theta_t), -sin(theta_t); sin(theta_t), cos(theta_t)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
                     'YData', rotated_satellite(2, :) + position(2));
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    p2 = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];
    set(h_link1, 'XData', [position(1), p1(1)], 'YData', [position(2), p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i)); % Update trajectory

    pause(0.05);
end
hold off;

% Extract components from state_history
heading_error = rad2deg(wrapToPi(state_history(3, :))); % Heading error in degrees
velocity_magnitude = vecnorm(state_history(4:5, :), 2, 1); % Velocity magnitude (m/s)

% Plot the results
figure;

% Plot Heading Error Over Time
subplot(3, 1, 1);
plot(time(1:size(state_history, 2)), heading_error, 'r', 'LineWidth', 1.5);
title('Heading Error Over Time');
xlabel('Time (s)');
ylabel('Heading Error (deg)');
grid on;

% Plot Velocity Over Time
subplot(3, 1, 2);
plot(time(1:size(state_history, 2)), velocity_magnitude, 'b', 'LineWidth', 1.5);
title('Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

% Plot Path Error Over Time
subplot(3, 1, 3);
plot(time(1:length(path_error_history)), path_error_history, 'k', 'LineWidth', 1.5);
title('Path Error Over Time');
xlabel('Time (s)');
ylabel('Path Error (m)');
grid on;


%%
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
gamma = 5;                % Adaptation gain for parameter updates
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target
velocity_limit = 1.0;     % Maximum translational velocity
omega_limit = 1.0;        % Maximum angular velocity
dq_limit = 5;             % Maximum joint velocity
q_limit = pi;             % Joint angle limit
damping_factor = 0.1;     % Damping factor for velocity and angular velocity


% Thruster force magnitudes
forward_thrust_magnitude = 0.1; % Maximum forward thrust

% Gains
k_theta = 0.2;    % Gain for heading error correction
k_l = 0.5;        % Gain for lateral deviation correction
k_v = 1;          % Gain for velocity control

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_t, v, omega, q1, q2, dq1, dq2]
state = [0; 0; pi/4; 0.1; 0; deg2rad(30); deg2rad(-45); 0; 0];

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length);
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);
heading_error_history = zeros(1, num_steps);
velocity_history = zeros(1, num_steps);
path_error_history = zeros(1, num_steps);

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Manipulator Dynamics');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Define spacecraft and manipulator shapes
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]';
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5); % Path tracking line

% Main simulation loop
for i = 1:num_steps
    % Extract state variables
    s = state(1); l = state(2); theta_t = state(3);
    v = state(4); omega = state(5);
    q1 = state(6); q2 = state(7); dq1 = state(8); dq2 = state(9);

    % Scale s to index the path array
    index = 1 + (s / path_length) * (path_resolution - 1);
    index = max(1, min(round(index), path_resolution)); % Clamp to valid range

    % Retrieve current position and tangent
    current_position = path(:, index);
    tangent = path_tangent(:, index);

    % Compute the satellite's position in global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;

    % Compute the heading error
    desired_heading = atan2(tangent(2), tangent(1));
    heading_error = wrapToPi(desired_heading - theta_t);
    heading_error_history(i) = rad2deg(heading_error); % Store heading error in degrees

    % Store velocity
    velocity_history(i) = v;

    % Compute path error (lateral deviation)
    path_error = l; % l is already the lateral deviation in Frenet-Serret dynamics
    path_error_history(i) = path_error;

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Control velocity and lateral deviation
    torque = -k_theta * theta_t; % Control heading error

    % Compute manipulator reaction torques
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];
    H = [m1_hat * l1^2 + m2_hat * (l1^2 + l2^2), m2_hat * l1 * l2 * cos(q2);
         m2_hat * l1 * l2 * cos(q2), m2_hat * l2^2];
    C = [0, -m2_hat * l1 * l2 * sin(q2) * dq2;
         m2_hat * l1 * l2 * sin(q2) * dq1, 0];
    tau_manipulator = -J' * (Kp_end_effector * [l; 0] + Kd_end_effector * [v; 0]); % End-effector torque

    ddq = pinv(H) * (tau_manipulator - C * [dq1; dq2]); % Joint accelerations

    % Enforce limits on joint velocities and angles
    dq1 = max(min(dq1 + ddq(1) * dt, dq_limit), -dq_limit);
    dq2 = max(min(dq2 + ddq(2) * dt, dq_limit), -dq_limit);
    q1 = max(min(q1 + dq1 * dt, q_limit), -q_limit);
    q2 = max(min(q2 + dq2 * dt, q_limit), -q_limit);

    % Apply manipulator torque to satellite angular dynamics
    tau_reaction = sum(tau_manipulator); % Reaction torque from manipulator
    angular_acceleration = (torque + tau_reaction) / inertia_satellite;

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Control velocity and lateral deviation
    torque = -k_theta * theta_t; % Control heading error

    % Compute manipulator reaction torques
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];
    H = [m1_hat * l1^2 + m2_hat * (l1^2 + l2^2), m2_hat * l1 * l2 * cos(q2);
         m2_hat * l1 * l2 * cos(q2), m2_hat * l2^2];
    C = [0, -m2_hat * l1 * l2 * sin(q2) * dq2;
         m2_hat * l1 * l2 * sin(q2) * dq1, 0];
    tau_manipulator = -J' * (Kp_end_effector * [l; 0] + Kd_end_effector * [v; 0]); % End-effector torque

    ddq = pinv(H) * (tau_manipulator - C * [dq1; dq2]); % Joint accelerations

    % Update joint states
    dq1 = dq1 + ddq(1) * dt;
    dq2 = dq2 + ddq(2) * dt;
    q1 = q1 + dq1 * dt;
    q2 = q2 + dq2 * dt;

    % Apply manipulator torque to satellite angular dynamics
    tau_reaction = -sum(tau_manipulator); % Reaction torque from manipulator
    angular_acceleration = (torque + tau_reaction) / inertia_satellite;

    % Update Frenet-Serret dynamics
    s_dot = v * cos(theta_t);
    l_dot = v * sin(theta_t);
    theta_t_dot = omega;
    v_dot = thrust / mass_satellite - damping_factor * v;
    omega_dot = angular_acceleration - damping_factor * omega;

    % Update the state
    state = state + [s_dot; l_dot; theta_t_dot; v_dot; omega_dot; dq1; dq2; ddq] * dt;

    % Update visualization
    R = [cos(theta_t), -sin(theta_t); sin(theta_t), cos(theta_t)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
                     'YData', rotated_satellite(2, :) + position(2));
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    p2 = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];
    set(h_link1, 'XData', [position(1), p1(1)], 'YData', [position(2), p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05);
end
hold off;

heading_error_history(i) = rad2deg(wrapToPi(theta_t));
velocity_history(i) = v;

% Plot Results
figure;

% Plot Heading Error Over Time
subplot(3, 1, 1);
plot(time(1:length(heading_error_history)), heading_error_history, 'r', 'LineWidth', 1.5);
title('Heading Error Over Time');
xlabel('Time (s)');
ylabel('Heading Error (deg)');
grid

% Plot Velocity Over Time
subplot(3, 1, 2);
plot(time(1:size(state_history, 2)), velocity_magnitude, 'b', 'LineWidth', 1.5);
title('Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;


% Plot Path Error Over Time
subplot(3, 1, 3);
plot(time(1:length(path_error_history)), path_error_history, 'k', 'LineWidth', 1.5);
title('Path Error Over Time');
xlabel('Time (s)');
ylabel('Path Error (m)');
grid on;




