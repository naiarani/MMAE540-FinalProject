%% Path following (min linear distance to target)

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
boundary_limit = 10;      % Bounds for the simulation space
target_position = [4; 4]; % Target position
target_velocity = 0.5;    % Desired velocity at the target
goal_tolerance = 0.1;     % Position tolerance for reaching the target

% Gains
k_theta = 2;    % Gain for heading error correction
k_v = 5;        % Gain for velocity correction
k_l = 1;        % Gain for lateral deviation correction
k_s = 2;        % Gain for progress correction

% Simulation parameters
t_final = 100; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_tilde, v, omega]
state = [0; 0; pi/4; 0.1; 0]; % Initial progress, lateral deviation, heading error, velocity, angular velocity

% Path definition (straight-line path to target)
path_start = [0; 0];
path_end = target_position;
path_direction = (path_end - path_start) / norm(path_end - path_start);

% Data storage for visualization
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);

% Define the satellite shape
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square satellite

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot([path_start(1), path_end(1)], [path_start(2), path_end(2)], 'k--', 'LineWidth', 1.5);

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Initialize satellite plot
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');

% Initialize trajectory plot
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Current states
    s = state(1); % Progress along path
    l = state(2); % Lateral deviation
    theta_tilda = state(3); % Heading error
    v = state(4); % Velocity
    w = state(5); % Angular velocity

    % Compute the current position of the spacecraft in global frame
    position = path_start + s * path_direction + [-path_direction(2); path_direction(1)] * l;
    position_history(:, i) = position;

    % Check if the spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance && abs(v - target_velocity) < 0.05
        disp('Spacecraft reached the target with desired velocity!');
        state_history = state_history(:, 1:i); % Trim unused steps
        position_history = position_history(:, 1:i); % Trim unused steps
        break;
    end

    % Control laws for Frenet-Serret dynamics
    w = -k_theta * theta_tilda; % Angular velocity to reduce heading error
    thrust = -k_l * l - k_s * (v - target_velocity); % Thrust to reduce lateral deviation and velocity error

    % Update Frenet-Serret dynamics
    theta_tilde_dot = w; % Change in heading error
    l_dot = v * sin(theta_tilda); % Change in lateral deviation
    s_dot = v * cos(theta_tilda); % Change in progress
    v_dot = thrust / mass_satellite; % Change in velocity

    % Update the state
    state = state + [s_dot; l_dot; theta_tilde_dot; v_dot; w] * dt;

    % Update animation
    R = [cos(theta_tilda), -sin(theta_tilda); sin(theta_tilda), cos(theta_tilda)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
        'YData', rotated_satellite(2, :) + position(2));
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05); % Animation speed
end
hold off;

% Additional plots
figure;

% Heading error over time
subplot(3, 1, 1);
plot(time(1:size(state_history, 2)), rad2deg(state_history(3, :)), 'r', 'LineWidth', 1.5);
title('Heading Error Over Time');
xlabel('Time (s)');
ylabel('Heading Error (deg)');
grid on;

% Velocity over time
subplot(3, 1, 2);
plot(time(1:size(state_history, 2)), state_history(4, :), 'b', 'LineWidth', 1.5);
title('Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

% Lateral deviation over time
subplot(3, 1, 3);
plot(time(1:size(state_history, 2)), state_history(2, :), 'k', 'LineWidth', 1.5);
title('Lateral Deviation Over Time');
xlabel('Time (s)');
ylabel('Lateral Deviation (m)');
grid on;

%% Sinuosoidal Path Following
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
target_velocity = 0.5;    % Desired velocity at the target
goal_tolerance = 0.1;     % Position tolerance for reaching the target

% Thruster force magnitudes
forward_thrust_magnitude = 0.2; % Maximum forward thrust
stopping_thrust_magnitude = 0.2; % Maximum stopping (opposing) thrust

% Gains
k_theta = 5;    % Gain for heading error correction
k_l = 2;        % Gain for lateral deviation correction
k_v = 10;       % Gain for velocity control

% Simulation parameters
t_final = 100; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [x, y, theta, vx, vy, omega]
state = [0; 0; 0; 0; 0; 0]; % Start at (0, 0) with theta = 0, no velocity or angular velocity

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length); % Sinusoidal shape
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage for visualization
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);
path_error_history = zeros(1, num_steps); % Store path error over time

% Define the satellite shape (square, initially at theta = 0)
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]';

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Orientation and Thrust Control');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Initialize satellite plot
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');

% Initialize trajectory plot
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Current states
    position = state(1:2); % Current position
    theta = state(3); % Current orientation
    velocity = state(4:5); % Current velocity
    omega = state(6); % Current angular velocity

    % Find the closest point on the path
    [~, closest_idx] = min(vecnorm(path - position));
    path_point = path(:, closest_idx); % Closest point on the path
    path_tan = path_tangent(:, closest_idx); % Path tangent at the closest point

    % Compute heading error: Angle between current orientation and path tangent
    desired_theta = atan2(path_tan(2), path_tan(1)); % Desired orientation along the path
    theta_error = wrapToPi(desired_theta - theta);

    % Compute lateral deviation: Perpendicular distance to the path
    path_error_vector = path_point - position; % Vector to the path point
    path_error = dot([-path_tan(2); path_tan(1)], path_error_vector); % Signed lateral deviation
    path_error_history(i) = path_error; % Store path error

    % Control law for angular velocity
    omega = k_theta * theta_error; % Proportional control for heading error

    % Control thrust
    % Forward thrust to correct velocity error and follow the path
    forward_thrust = k_v * (target_velocity - norm(velocity)) * path_tan;

    % Lateral correction thrust to minimize path error
    lateral_correction = -k_l * path_error * [-path_tan(2); path_tan(1)];

    % Total thrust
    total_thrust = forward_thrust + lateral_correction;

    % Saturate thrust to maximum values
    thrust_magnitude = norm(total_thrust);
    if thrust_magnitude > forward_thrust_magnitude
        total_thrust = (total_thrust / thrust_magnitude) * forward_thrust_magnitude;
    end

    acceleration = total_thrust / mass_satellite;

    % Update dynamics
    velocity = velocity + acceleration * dt;
    position = position + velocity * dt;

    angular_acceleration = omega / inertia_satellite;
    omega = omega + angular_acceleration * dt;
    theta = theta + omega * dt;

    % Update the state
    state = [position; theta; velocity; omega];

    % Check if the spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance && abs(norm(velocity) - target_velocity) < 0.05
        disp('Spacecraft reached the target with desired velocity and orientation!');
        state_history = state_history(:, 1:i); % Trim unused steps
        position_history = position_history(:, 1:i); % Trim unused steps
        path_error_history = path_error_history(1:i); % Trim unused steps
        break;
    end

    % Update animation
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
        'YData', rotated_satellite(2, :) + position(2));
    position_history(:, i) = position; % Store position history
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05); % Animation speed
end
hold off;

% Additional plots
figure;

% Path error over time
subplot(3, 1, 1);
plot(time(1:length(path_error_history)), path_error_history, 'b', 'LineWidth', 1.5);
title('Path Error Over Time');
xlabel('Time (s)');
ylabel('Path Error (m)');
grid on;

% Heading error over time
subplot(3, 1, 2);
plot(time(1:size(state_history, 2)), rad2deg(wrapToPi(state_history(3, :))), 'r', 'LineWidth', 1.5);
title('Heading Error Over Time');
xlabel('Time (s)');
ylabel('Heading Error (deg)');
grid on;

% Velocity magnitude over time
subplot(3, 1, 3);
plot(time(1:size(state_history, 2)), vecnorm(state_history(4:5, :), 2, 1), 'g', 'LineWidth', 1.5);
title('Velocity Magnitude Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;


%% Unstable sine wave path following

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target
kappa = 2 * pi / 10;      % Curvature for sinusoidal path (wavelength = 10)

% Gains
k_s = 2;    % Progress control gain
k_l = 10;    % Lateral deviation control gain
k_theta = 0.5; % Heading control gain
k_v = 1;    % Velocity control gain

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [s, l, theta_t, v, omega]
state = [0; 0; pi/4; 0.1; 0]; % Initial progress, lateral deviation, heading error, velocity, angular velocity

% Path definition (sinusoidal path)
path_length = 10; % Path length
path = @(s) [s; 2 * sin(kappa * s)]; % Sinusoidal path as a function of progress (s)

% Data storage for visualization
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);
path_error_history = zeros(1, num_steps);

% Define the satellite shape
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square satellite

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Thrust Boosts and Orientation Control');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
s_path = linspace(0, path_length, 100);
path_points = arrayfun(path, s_path, 'UniformOutput', false);
path_points = cat(2, path_points{:});
plot(path_points(1, :), path_points(2, :), 'k--', 'LineWidth', 1.5);

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Initialize satellite plot
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');

% Initialize trajectory plot
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Extract Frenet-Serret states
    s = state(1); % Progress along the path
    l = state(2); % Lateral deviation
    theta_t = state(3); % Heading error
    v = state(4); % Velocity
    omega = state(5); % Angular velocity

    % Compute the current position and tangent on the path
    current_position = path(s);
    tangent = [1; 2 * kappa * cos(kappa * s)];
    tangent = tangent / norm(tangent); % Normalize tangent

    % Compute the spacecraft's position in the global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;
    
    % Calculate the path error (lateral deviation)
    path_error_vector = position - current_position; % Vector from path point to spacecraft
    path_error = dot([-tangent(2); tangent(1)], path_error_vector); % Signed lateral deviation
    path_error_history(i) = path_error; % Store the path error

    % Compute the global position of the spacecraft
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;

    % Check if the spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance && abs(v - target_velocity) < 0.05
        disp('Spacecraft reached the target with desired velocity!');
        state_history = state_history(:, 1:i); % Trim unused steps
        position_history = position_history(:, 1:i); % Trim unused steps
        break;
    end

    % Control laws for Frenet-Serret dynamics
    thrust = k_s * (target_velocity - v) - k_l * l; % Thrust to control velocity and lateral deviation
    torque = k_theta * theta_t; % Torque to correct heading error

    % Apply thrust bursts
    if abs(l) > 0.05 || abs(theta_t) > deg2rad(5)
        acceleration = thrust / mass_satellite;
    else
        acceleration = 0; % Stop thrusting when aligned
    end

    angular_acceleration = torque / inertia_satellite;

    % Update Frenet-Serret dynamics
    theta_t_dot = omega; % Change in heading error
    l_dot = v * sin(theta_t); % Change in lateral deviation
    s_dot = v * cos(theta_t); % Change in progress
    v_dot = acceleration; % Change in velocity
    omega_dot = angular_acceleration; % Change in angular velocity

    % Update the state
    state = state + [s_dot; l_dot; theta_t_dot; v_dot; omega_dot] * dt;

    % Update animation
    R = [cos(state(3)), -sin(state(3)); sin(state(3)), cos(state(3))];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
        'YData', rotated_satellite(2, :) + position(2));
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05); % Animation speed
end
hold off;

% Controllability Analysis
disp('Checking Controllability...');
A = [0 1 0 0 0; 
     0 0 -1 0 0;
     0 0 0 1 0;
     0 0 0 0 1;
     0 0 0 0 -k_theta]; % Linearized A matrix (example structure)
B = [0; 0; 0; 1/mass_satellite; 1/inertia_satellite]; % Linearized B matrix
controllability_matrix = ctrb(A, B);
rank_controllability = rank(controllability_matrix);
disp(['Controllability Matrix Rank: ', num2str(rank_controllability)]);
disp(['Number of States: ', num2str(size(A, 1))]);

if rank_controllability == size(A, 1)
    disp('The system is controllable.');
else
    disp('The system is NOT controllable.');
end

% Extract components from `state_history`
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


%% Fernet Serret Path Following ! (discritized path)
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
goal_tolerance = 0.1;     % Position tolerance for reaching the target
target_velocity = 0.5;    % Desired velocity at the target

% Thruster force magnitudes
forward_thrust_magnitude = 0.1; % Maximum forward thrust
stopping_thrust_magnitude = 0.1; % Maximum stopping (opposing) thrust

% Gains
k_theta = 0.2;    % Gain for heading error correction
k_l = 0.5;        % Gain for lateral deviation correction
k_v = 1;       % Gain for velocity control

% Simulation parameters
t_final = 50; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial Frenet-Serret state: [s, l, theta_t, v, omega]
state = [0; 0; pi/4; 0.1; 0]; % Initial progress, lateral deviation, heading error, velocity, angular velocity

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

% Data storage for visualization
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);
path_error_history = zeros(1, num_steps);

% Define the satellite shape
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square satellite

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Frenet-Serret States');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Initialize satellite plot
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');

% Initialize trajectory plot
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Extract Frenet-Serret states
    s = state(1); % Progress along the path
    l = state(2); % Lateral deviation
    theta_t = state(3); % Heading error
    v = state(4); % Velocity
    omega = state(5); % Angular velocity

    % Scale s to index the path array
    index = 1 + (s / path_length) * (path_resolution - 1);
    index = max(1, min(round(index), path_resolution)); % Clamp to valid range  

    % Retrieve the current position and tangent on the path
    current_position = path(:, index);
    tangent = path_tangent(:, index);
    
    % Compute the spacecraft's position in the global frame
    position = current_position + [-tangent(2); tangent(1)] * l; % Offset from path
    position_history(:, i) = position;
    
    % Calculate the path error (lateral deviation)
    path_error_vector = position - current_position; % Vector from path point to spacecraft
    path_error = dot([-tangent(2); tangent(1)], path_error_vector); % Signed lateral deviation
    path_error_history(i) = path_error; % Store the path error

    % Check if the spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance && abs(v - target_velocity) < 0.05
        disp('Spacecraft reached the target with desired velocity and orientation!');
        state_history = state_history(:, 1:i); % Trim unused steps
        position_history = position_history(:, 1:i); % Trim unused steps
        break;
    end

    % Control laws for Frenet-Serret dynamics
    thrust = k_v * (target_velocity - v) - k_l * l; % Thrust to control velocity and lateral deviation
    torque = -k_theta * theta_t; % Torque to correct heading error

    % Apply thrust bursts
    if abs(l) > 0.05 || abs(theta_t) > deg2rad(5)
        acceleration = thrust / mass_satellite;
    else
        acceleration = 0; % Stop thrusting when aligned
    end

    angular_acceleration = torque / inertia_satellite;

    % Update Frenet-Serret dynamics
    theta_t_dot = omega; % Change in heading error
    l_dot = v * sin(theta_t); % Change in lateral deviation
    s_dot = v * cos(theta_t); % Change in progress
    v_dot = acceleration; % Change in velocity
    omega_dot = angular_acceleration; % Change in angular velocity

    % Update the state
    state = state + [s_dot; l_dot; theta_t_dot; v_dot; omega_dot] * dt;

    % Update animation
    R = [cos(state(3)), -sin(state(3)); sin(state(3)), cos(state(3))];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
        'YData', rotated_satellite(2, :) + position(2));
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05); % Animation speed
end
hold off;

% Ensure your simulation results are already calculated in:
% `state_history` (state trajectory),
% `time` (time vector), and
% `path_error_history` (path error over time).

% Extract components from `state_history`
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

%% Semi-decent path following but with global frame instead of fernet 
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
boundary_limit = 10;      % Bounds for the simulation space
target_position = [10; 0]; % Target position
target_velocity = 0.5;    % Desired velocity at the target
goal_tolerance = 0.1;     % Position tolerance for reaching the target

% Thruster force magnitudes
forward_thrust_magnitude = 0.2; % Maximum forward thrust
stopping_thrust_magnitude = 0.2; % Maximum stopping (opposing) thrust

% Gains
k_theta = 5;    % Gain for heading error correction
k_l = 2;        % Gain for lateral deviation correction
k_v = 10;       % Gain for velocity control

% Simulation parameters
t_final = 100; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial state: [x, y, theta, vx, vy, omega]
state = [0; 0; pi/4; 0.1; 0; 0]; % Initial position, orientation, velocity, and angular velocity

% Path definition: sinusoidal path to the target
path_length = target_position(1);
path_resolution = 100;
path_x = linspace(0, path_length, path_resolution);
path_y = 2 * sin(2 * pi * path_x / path_length); % Sinusoidal shape
path = [path_x; path_y];

% Compute tangent vectors for the path
path_tangent_x = gradient(path(1, :));
path_tangent_y = gradient(path(2, :));
path_tangent = [path_tangent_x; path_tangent_y];
path_tangent = path_tangent ./ vecnorm(path_tangent); % Normalize tangents

% Data storage for visualization
state_history = zeros(length(state), num_steps);
position_history = zeros(2, num_steps);
path_error_history = zeros(1, num_steps); % Store path error over time

% Define the satellite shape
satellite_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square satellite

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft Path Following with Orientation Control');
xlabel('X Position');
ylabel('Y Position');

% Draw the path
plot(path(1, :), path(2, :), 'k--', 'LineWidth', 1.5);

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Initialize satellite plot
h_satellite = fill(satellite_shape(1, :), satellite_shape(2, :), 'b');

% Initialize trajectory plot
h_trajectory = plot(NaN, NaN, 'g', 'LineWidth', 1.5);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Current states
    position = state(1:2); % Current position
    theta = state(3); % Current orientation
    velocity = state(4:5); % Current velocity
    omega = state(6); % Current angular velocity

    % Find the closest point on the path
    [~, closest_idx] = min(vecnorm(path - position));
    path_point = path(:, closest_idx); % Closest point on the path
    path_tan = path_tangent(:, closest_idx); % Path tangent at the closest point

    % Compute heading error: Angle between current orientation and path tangent
    desired_theta = atan2(path_tan(2), path_tan(1)); % Desired orientation along the path
    theta_error = wrapToPi(desired_theta - theta);

    % Compute lateral deviation: Perpendicular distance to the path
    path_error_vector = path_point - position; % Vector to the path point
    path_error = dot([-path_tan(2); path_tan(1)], path_error_vector); % Signed lateral deviation
    path_error_history(i) = path_error; % Store path error

    % Control law for angular velocity
    omega = k_theta * theta_error; % Proportional control for heading error

    % Control thrust: Apply forward thrust if aligned with the path
    if abs(theta_error) < deg2rad(10) % Only thrust if aligned within a threshold
        forward_thrust = forward_thrust_magnitude * [cos(theta); sin(theta)];
    else
        forward_thrust = [0; 0]; % No forward thrust if not aligned
    end

    % Apply stopping thrust if velocity exceeds the target or lateral deviation is too high
    stopping_thrust = [0; 0];
    if norm(velocity) > target_velocity || abs(path_error) > goal_tolerance
        stopping_thrust = -stopping_thrust_magnitude * (velocity / norm(velocity)); % Opposes current velocity
    end

    % Total thrust
    total_thrust = forward_thrust + stopping_thrust;
    acceleration = total_thrust / mass_satellite;

    % Update dynamics
    velocity = velocity + acceleration * dt;
    position = position + velocity * dt;

    angular_acceleration = omega / inertia_satellite;
    omega = omega + angular_acceleration * dt;
    theta = theta + omega * dt;

    % Update the state
    state = [position; theta; velocity; omega];

    % Check if the spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance && abs(norm(velocity) - target_velocity) < 0.05
        disp('Spacecraft reached the target with desired velocity and orientation!');
        state_history = state_history(:, 1:i); % Trim unused steps
        position_history = position_history(:, 1:i); % Trim unused steps
        path_error_history = path_error_history(1:i); % Trim unused steps
        break;
    end

    % Update animation
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated_satellite = R * satellite_shape;
    set(h_satellite, 'XData', rotated_satellite(1, :) + position(1), ...
        'YData', rotated_satellite(2, :) + position(2));
    position_history(:, i) = position; % Store position history
    set(h_trajectory, 'XData', position_history(1, 1:i), 'YData', position_history(2, 1:i));

    pause(0.05); % Animation speed
end
hold off;

% Additional plots
figure;

% Heading error over time
subplot(3, 1, 1);
plot(time(1:size(state_history, 2)), rad2deg(wrapToPi(state_history(3, :) - desired_theta)), 'r', 'LineWidth', 1.5);
title('Heading Error Over Time');
xlabel('Time (s)');
ylabel('Heading Error (deg)');
grid on;

% Velocity over time
subplot(3, 1, 2);
plot(time(1:size(state_history, 2)), vecnorm(state_history(4:5, :), 2, 1), 'b', 'LineWidth', 1.5);
title('Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

% Path error over time
subplot(3, 1, 3);
plot(time(1:length(path_error_history)), path_error_history, 'k', 'LineWidth', 1.5);
title('Path Error Over Time');
xlabel('Time (s)');
ylabel('Path Error (m)');
grid on;
