% 540 Final Project

clc;
clear;
close all;

% Constants and Setup
mass_satellite = 10.0;    % Mass of the satellite for acceleration calculations
inertia_satellite = 10.0; % Moment of inertia of the satellite
target_position = [4; 4]; % Target position in 2D space
goal_tolerance = 0.1;     % Distance tolerance for reaching the goal

% Control Gains
Kp_orientation = 0.1;
Kd_orientation = 0.05;
Kp_thruster = 0.5;

% Simulation Parameters
t_final = 20;             % Total simulation time
dt = 0.1;                 % Time step
num_steps = t_final / dt;

% Initial State: [x, y, theta, vx, vy, omega]
state = [0; 0; 0; 0; 0; 0];  % Starting at the origin

% Array to store the state history
state_history = zeros(6, num_steps);

for i = 1:num_steps
    % Store the current state
    state_history(:, i) = state;

    % Extract current position and orientation
    position = state(1:2);
    theta = state(3);

    % Check if the spacecraft has reached the goal
    if norm(position - target_position) < goal_tolerance
        disp('Target reached!');
        state_history = state_history(:, 1:i); % Trim unused steps
        break;
    end

    % PD Control for CMG to adjust orientation
    position_error = target_position - position;
    control_input = Kp_orientation * norm(position_error) - Kd_orientation * theta;

    % Thruster control to move toward the target position
    thrust_force = Kp_thruster * position_error;

    % Update dynamics with thrust and control input
    dposition_dt = state(4:5);                       % Velocity affects position
    dtheta_dt = state(6);                            % Angular velocity affects orientation
    dvelocity_dt = thrust_force / mass_satellite;    % Thrusters affect velocity
    domega_dt = control_input / inertia_satellite;   % Control input affects angular velocity

    % Update the state based on dynamics
    state(1:2) = state(1:2) + dposition_dt * dt;     % Update position
    state(3) = state(3) + dtheta_dt * dt;            % Update orientation
    state(4:5) = state(4:5) + dvelocity_dt * dt;     % Update velocity
    state(6) = state(6) + domega_dt * dt;            % Update angular velocity
end

% Visualization
figure;
hold on;
axis equal;
xlim([-5, 10]);
ylim([-5, 10]);
title('Spacecraft Moving Towards Target');
xlabel('X Position');
ylabel('Y Position');

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Define the spacecraft as a triangular polygon for rotation
spacecraft_shape = [0.25, 0; -0.25, 0.15; -0.25, -0.15]'; % Shape vertices
h_spacecraft = fill(spacecraft_shape(1, :), spacecraft_shape(2, :), 'b');

for k = 1:size(state_history, 2)
    % Extract the current position and orientation
    x = state_history(1, k);
    y = state_history(2, k);
    theta = state_history(3, k);
    
    % Rotation matrix for the orientation angle
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated_shape = R * spacecraft_shape; % Rotate the spacecraft shape

    % Update the spacecraft position
    set(h_spacecraft, 'XData', rotated_shape(1, :) + x, 'YData', rotated_shape(2, :) + y);

    pause(0.05); % Control the animation speed
end
hold off;


% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Draw the spacecraft as a rectangle and update its position in each step
spacecraft_shape = rectangle('Position', [0, 0, 0.5, 0.3], 'FaceColor', 'b');

for k = 1:size(state_history, 2)
    % Update the spacecraft position and orientation
    x = state_history(1, k);
    y = state_history(2, k);
    theta = state_history(3, k);
    
    % Set spacecraft position and rotation
    spacecraft_shape.Position = [x - 0.25, y - 0.15, 0.5, 0.3];
    spacecraft_shape.Rotation = rad2deg(theta); % Rotate to orientation angle

    pause(0.05); % Control the animation speed
end
hold off;

%%

clc;
clear;
close all;

% Constants for spacecraft and manipulator
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
m1 = 5; m2 = 5;           % Link masses
lc1 = l1 / 2; lc2 = l2 / 2; % Link center of mass
I1 = (1/12) * m1 * l1^2;  % Inertia of link 1
I2 = (1/12) * m2 * l2^2;  % Inertia of link 2

% Control Gains for manipulator and spacecraft
Kp_orientation = 0.1; Kd_orientation = 0.05;
Kp_manipulator = 60; Kd_manipulator = 20;
Kp_thruster = 20;

% Simulation parameters
t_final = 20; dt = 0.1;
num_steps = t_final / dt;

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Desired joint angles for manipulator
q1_desired = deg2rad(60); q2_desired = deg2rad(-30);
dq1_desired = 0; dq2_desired = 0;

% Store state history for visualization
state_history = zeros(length(state), num_steps);

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

for i = 1:num_steps
    % Store current state
    state_history(:, i) = state;
    
    % Extract current values
    position = state(1:2);
    theta = state(3);
    q1 = state(4); q2 = state(5);
    velocity = state(6:7);
    omega = state(8);
    dq1 = state(9); dq2 = state(10);
    
    % Check if spacecraft has reached the target
    if norm(position - target_position) < goal_tolerance
        disp('Target reached!');
        state_history = state_history(:, 1:i); % Trim unused steps
        break;
    end
    
    % Spacecraft PD control for orientation
    position_error = target_position - position;
    tau_spacecraft = Kp_orientation * norm(position_error) - Kd_orientation * theta;

    % Manipulator PD control for joint angles
    tau1 = Kp_manipulator * (q1_desired - q1) + Kd_manipulator * (dq1_desired - dq1);
    tau2 = Kp_manipulator * (q2_desired - q2) + Kd_manipulator * (dq2_desired - dq2);
    tau_manipulator = [tau1; tau2];
    
    % Dynamics of manipulator
    H11 = m1*lc1^2 + I1 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(q2)) + I2;
    H22 = m2*lc2^2 + I2;
    H12 = m2*l1*lc2*cos(q2) + I2;
    H = [H11, H12; H12, H22];
    
    h = m2 * l1 * lc2 * sin(q2);
    C = [-h*dq2, -h*(dq1+dq2); h*dq1, 0];
    
    % Compute angular accelerations for the manipulator
    ddq = H \ (tau_manipulator - C * [dq1; dq2]);
    
    % Update manipulator joint velocities and positions
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Compute spacecraft dynamics
    thrust_force = Kp_thruster * position_error; % Thruster control
    dposition_dt = velocity; % Velocity affects position
    dvelocity_dt = thrust_force / mass_satellite; % Force affects velocity
    domega_dt = tau_spacecraft / inertia_satellite; % Torque affects angular velocity

    % Update spacecraft state
    velocity = velocity + dvelocity_dt * dt;
    position = position + dposition_dt * dt;
    omega = omega + domega_dt * dt;
    theta = theta + omega * dt;

    % Update state vector
    state = [position; theta; q1; q2; velocity; omega; dq1; dq2];
end

% Visualization
figure;
hold on;
axis equal;
xlim([-5, 10]);
ylim([-5, 10]);
title('Spacecraft and 2-Link Manipulator');
xlabel('X Position');
ylabel('Y Position');

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Define spacecraft and manipulator shapes
spacecraft_shape = [0.25, 0; -0.25, 0.15; -0.25, -0.15]'; % Triangular spacecraft
link1 = [0, 0; l1, 0]'; % Link 1
link2 = [0, 0; l2, 0]'; % Link 2

% Create plots
h_spacecraft = fill(spacecraft_shape(1, :), spacecraft_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 2

for k = 1:size(state_history, 2)
    % Extract states
    x = state_history(1, k);
    y = state_history(2, k);
    theta = state_history(3, k);
    q1 = state_history(4, k);
    q2 = state_history(5, k);
    
    % Rotate and translate spacecraft
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated_spacecraft = R * spacecraft_shape;
    set(h_spacecraft, 'XData', rotated_spacecraft(1, :) + x, 'YData', rotated_spacecraft(2, :) + y);
    
    % Forward kinematics for manipulator
    p1 = [x; y] + R * [l1*cos(q1); l1*sin(q1)];
    p2 = p1 + R * [l2*cos(q1 + q2); l2*sin(q1 + q2)];
    set(h_link1, 'XData', [x, p1(1)], 'YData', [y, p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);

    pause(0.05); % Control animation speed
end
hold off;

%%
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
l1 = 1; l2 = 1;           % Link lengths
m1 = 5; m2 = 5;           % Link masses
Kp_end_effector = 50; Kd_end_effector = 10; % End-effector PD gains
thruster_force = 0.5;      % Small thrust boost
boundary_limit = 10;       % Bounds for the simulation space

% Simulation parameters
t_final = 200; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

% Data storage for analysis and visualization
state_history = zeros(length(state), num_steps);
thrust_history = zeros(2, num_steps); % To track thruster activations
x_tilde_history = zeros(2, num_steps); % End-effector error

for i = 1:num_steps
    % Store state
    state_history(:, i) = state;

    % Extract current values
    position = state(1:2);
    theta = state(3);
    q1 = state(4); q2 = state(5);
    velocity = state(6:7);
    omega = state(8);
    dq1 = state(9); dq2 = state(10);

    % Compute current end-effector position using forward kinematics
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    end_effector = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];

    % Check if the end effector has reached the target
    if norm(end_effector - target_position) < goal_tolerance
        disp('End effector reached the target!');
        state_history = state_history(:, 1:i); % Trim unused steps
        thrust_history = thrust_history(:, 1:i); % Trim thrust history
        x_tilde_history = x_tilde_history(:, 1:i); % Trim error history
        break;
    end

    % End-effector control
    x_tilde = end_effector - target_position; % Position error
    dx_tilde = R * ([l1 * -sin(q1), -l2 * sin(q1 + q2); l1 * cos(q1), l2 * cos(q1 + q2)] * [dq1; dq2]); % Velocity error
    x_tilde_history(:, i) = x_tilde; % Log position error

    % Jacobian matrix for the 2-link manipulator
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];

    % PD control for the manipulator in end-effector space
    tau_manipulator = -J' * (Kp_end_effector * x_tilde + Kd_end_effector * dx_tilde);

    % Dynamics of manipulator
    ddq = [0; 0]; % Manipulator dynamics placeholder (or compute as needed)

    % Update manipulator state
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Thruster control: Constant thrust toward target
    position_error = target_position - position; % Error in spacecraft position
    thrust = thruster_force * position_error / norm(position_error); % Normalize thrust vector
    thrust_history(:, i) = thrust; % Log thrust

    % Compute spacecraft dynamics
    dposition_dt = velocity; % Velocity affects position
    dvelocity_dt = thrust / mass_satellite; % Thrusters affect velocity
    position = position + dposition_dt * dt;
    velocity = velocity + dvelocity_dt * dt;

    % Update state vector
    state = [position; theta; q1; q2; velocity; omega; dq1; dq2];
end

% Visualization
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft and 2-Link Manipulator with Thrusters');
xlabel('X Position');
ylabel('Y Position');

% Draw target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Define spacecraft and manipulator shapes
spacecraft_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square spacecraft
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_spacecraft = fill(spacecraft_shape(1, :), spacecraft_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 2

for k = 1:size(state_history, 2)
    % Extract states
    x = state_history(1, k);
    y = state_history(2, k);
    theta = state_history(3, k);
    q1 = state_history(4, k);
    q2 = state_history(5, k);

    % Rotate and translate spacecraft
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated_spacecraft = R * spacecraft_shape;
    set(h_spacecraft, 'XData', rotated_spacecraft(1, :) + x, 'YData', rotated_spacecraft(2, :) + y);

    % Forward kinematics for manipulator
    p1 = [x; y] + R * [l1*cos(q1); l1*sin(q1)];
    p2 = p1 + R * [l2*cos(q1 + q2); l2*sin(q1 + q2)];
    set(h_link1, 'XData', [x, p1(1)], 'YData', [y, p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);

    pause(0.05); % Control animation speed
end
hold off;

% Additional plots
figure;
plot(time, vecnorm(x_tilde_history, 2, 1));
title('End-Effector Error Over Time');
xlabel('Time (s)');
ylabel('Error (m)');
