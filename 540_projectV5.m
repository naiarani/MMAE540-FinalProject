%% Adaptive Controller Tracking a moving target

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
thruster_force = 5;     % Thruster boost
boundary_limit = 10;      % Bounds for the simulation space

% Limits
q_limit = deg2rad(150);    % Maximum joint angle limit
dq_limit = deg2rad(50);    % Maximum joint velocity limit
omega_limit = 1;           % Maximum satellite angular velocity

% Simulation parameters
t_final = 150; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Target position and velocity
target_position = [4; 4]; % Initial position
target_velocity = [2; 1]; % Constant velocity (m/s)

goal_tolerance = 0.5;

% Data storage for analysis and visualization
state_history = zeros(length(state), num_steps);
x_tilde_history = zeros(2, num_steps); % End-effector error
m1_hat_history = zeros(1, num_steps); % Adaptive parameter m1_hat
m2_hat_history = zeros(1, num_steps); % Adaptive parameter m2_hat
target_position_history = zeros(2, num_steps); % Target position history

% Visualization setup
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft and 2-Link Manipulator with Moving Target');
xlabel('X Position');
ylabel('Y Position');

% Draw the target as a circle
h_target = plot(target_position(1), target_position(2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% Define spacecraft and manipulator shapes
spacecraft_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square spacecraft
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_spacecraft = fill(spacecraft_shape(1, :), spacecraft_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);

for i = 1:num_steps
    % Store state and target position
    state_history(:, i) = state;
    target_position_history(:, i) = target_position;
    m1_hat_history(i) = m1_hat;
    m2_hat_history(i) = m2_hat;

    % Extract current values
    position = state(1:2);
    phi = state(3); % Current orientation
    q1 = state(4); q2 = state(5);
    velocity = state(6:7);
    omega = state(8);
    dq1 = state(9); dq2 = state(10);

    % Compute current end-effector position using forward kinematics
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    end_effector = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];

    % Check if the end effector has reached the target
    if norm(end_effector - target_position) < goal_tolerance
        disp('End effector reached the target!');
        state_history = state_history(:, 1:i); % Trim unused steps
        x_tilde_history = x_tilde_history(:, 1:i); % Trim error history
        m1_hat_history = m1_hat_history(1:i);
        m2_hat_history = m2_hat_history(1:i);
        target_position_history = target_position_history(:, 1:i);
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

    % Adaptive control law
    m1_hat_dot = -gamma * (l1 * cos(q1)) * x_tilde(1);
    m2_hat_dot = -gamma * (l2 * cos(q1 + q2)) * x_tilde(1);

    % Update adaptive parameters
    m1_hat = max(m1_hat + m1_hat_dot * dt, 0); % Ensure non-negative
    m2_hat = max(m2_hat + m2_hat_dot * dt, 0); % Ensure non-negative

    % Reaction torque from manipulator affects spacecraft
    tau_reaction = -sum(tau_manipulator); % Reaction torque on spacecraft

    % Dynamics of manipulator
    H11 = m1_hat * 0.5^2 + m2_hat * (l1^2 + l2^2) + 2 * m2_hat * l1 * l2 * cos(q2) + inertia_satellite;
    H22 = m2_hat * l2^2 + inertia_satellite;
    H12 = m2_hat * l1 * l2 * cos(q2);
    H = [H11, H12; H12, H22];

    % Velocity-dependent forces (Coriolis)
    h = -m2_hat * l1 * l2 * sin(q2);
    C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];

    % Joint accelerations
    ddq = H \ (tau_manipulator - C * [dq1; dq2]);

    % Update manipulator state with limits
    dq1 = max(min(dq1 + ddq(1) * dt, dq_limit), -dq_limit);
    dq2 = max(min(dq2 + ddq(2) * dt, dq_limit), -dq_limit);
    q1 = max(min(q1 + dq1 * dt, q_limit), -q_limit);
    q2 = max(min(q2 + dq2 * dt, q_limit), -q_limit);

    % Update spacecraft dynamics using manipulator reaction torque
    domega_dt = tau_reaction / inertia_satellite; % Angular acceleration
    omega = max(min(omega + domega_dt * dt, omega_limit), -omega_limit); % Update angular velocity
    phi = phi + omega * dt; % Update orientation

    % Thruster control: Constant thrust toward target
    position_error = target_position - position; % Error in spacecraft position
    thrust = thruster_force * position_error / max(norm(position_error), 1e-6); % Normalize thrust vector

    % Compute spacecraft translational dynamics
    dposition_dt = velocity; % Velocity affects position
    dvelocity_dt = thrust / mass_satellite; % Thrusters affect velocity
    position = position + dposition_dt * dt;
    velocity = velocity + dvelocity_dt * dt;

    % Update state vector
    state = [position; phi; q1; q2; velocity; omega; dq1; dq2];

    % Update target position and ensure it stays within bounds
    target_position = target_position + target_velocity * dt;
    if any(abs(target_position) > boundary_limit)
        target_velocity = -target_velocity; % Reverse velocity if hitting a boundary
    end

    % Animation
    set(h_target, 'XData', target_position(1), 'YData', target_position(2));
    rotated_spacecraft = R * spacecraft_shape;
    set(h_spacecraft, 'XData', rotated_spacecraft(1, :) + position(1), 'YData', rotated_spacecraft(2, :) + position(2));
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    p2 = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];
    set(h_link1, 'XData', [position(1), p1(1)], 'YData', [position(2), p1(2)]);
    set(h_link2, 'XData', [p1(1), p2(1)], 'YData', [p1(2), p2(2)]);
    pause(0.05); % Adjust for animation speed
end

