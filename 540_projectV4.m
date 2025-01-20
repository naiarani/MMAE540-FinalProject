% MPC
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
boundary_limit = 10;      % Bounds for the simulation space
thruster_force = 0.5;

% MPC Parameters
horizon = 10; % Prediction horizon
dt = 0.1; % Time step
Q = eye(2) * 100; % State error weight
R = eye(2) * 1; % Control effort weight
thruster_penalty = 10; % Penalize excessive thruster use

% Limits
tau_limit = 5; % Maximum torque
q_limit = deg2rad(150);    % Maximum joint angle limit
dq_limit = deg2rad(50);    % Maximum joint velocity limit

% Simulation parameters
t_final = 100; 
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

% Data storage for analysis and visualization
state_history = zeros(length(state), num_steps);
x_tilde_history = zeros(2, num_steps); % End-effector error
tau_history = zeros(2, num_steps); % Torque history

for i = 1:num_steps
    % Store state
    state_history(:, i) = state;

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
        tau_history = tau_history(:, 1:i); % Trim torque history
        break;
    end
    
    % Compute end-effector error
    x_tilde = target_position - end_effector; % Error in end-effector position
    x_tilde_history(:, i) = x_tilde; % Store error history for plotting


    % Define the MPC optimization problem
    u_init = zeros(2 * horizon, 1); % Initial guess for torques
    lb = -tau_limit * ones(2 * horizon, 1); % Lower bound
    ub = tau_limit * ones(2 * horizon, 1); % Upper bound
    cost_function = @(u) mpc_cost_coupled(u, state, horizon, target_position, Q, R, thruster_penalty, l1, l2, dt);
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

    % Solve the MPC problem
    u_opt = fmincon(cost_function, u_init, [], [], [], [], lb, ub, [], options);

    % Extract the first control input
    tau = u_opt(1:2);
    tau_history(:, i) = tau; % Store torques

    % Dynamics of manipulator
    H11 = 5 * 0.5^2 + 5 * (l1^2 + l2^2) + 2 * 5 * l1 * l2 * cos(q2);
    H22 = 5 * l2^2;
    H12 = 5 * l1 * l2 * cos(q2);
    H = [H11, H12; H12, H22];
    h = -5 * l1 * l2 * sin(q2);
    C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];
    ddq = H \ (tau - C * [dq1; dq2]);

    % Reaction forces and torques on the spacecraft
    reaction_force = -R * [tau(1); tau(2)];
    domega_dt = sum(tau) / inertia_satellite; % Reaction torque
    omega = omega + domega_dt * dt;

    % Update manipulator state with limits
    dq1 = max(min(dq1 + ddq(1) * dt, dq_limit), -dq_limit);
    dq2 = max(min(dq2 + ddq(2) * dt, dq_limit), -dq_limit);
    q1 = max(min(q1 + dq1 * dt, q_limit), -q_limit);
    q2 = max(min(q2 + dq2 * dt, q_limit), -q_limit);

    % Thruster control: Thrusters are less used as reaction forces assist
    position_error = target_position - position; % Error in spacecraft position
    thrust = thruster_force * position_error / max(norm(position_error), 1e-6) - reaction_force;

    % Compute spacecraft translational dynamics
    dposition_dt = velocity; % Velocity affects position
    dvelocity_dt = thrust / mass_satellite; % Thrusters affect velocity
    position = position + dposition_dt * dt;
    velocity = velocity + dvelocity_dt * dt;

    % Update state vector
    state = [position; phi; q1; q2; velocity; omega; dq1; dq2];
end

% Visualization
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft and 2-Link Manipulator with Reaction Coupling');
xlabel('X Position');
ylabel('Y Position');

% Draw the target as a circle
viscircles(target_position', goal_tolerance, 'Color', 'r', 'LineWidth', 0.5);

% Define spacecraft and manipulator shapes
spacecraft_shape = [-0.5, -0.5; 0.5, -0.5; 0.5, 0.5; -0.5, 0.5]'; % Square spacecraft
link1 = [0, 0; l1, 0]';
link2 = [0, 0; l2, 0]';

% Create plots
h_spacecraft = fill(spacecraft_shape(1, :), spacecraft_shape(2, :), 'b');
h_link1 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 1
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2);

for k = 1:size(state_history, 2)
    % Extract states
    x = state_history(1, k);
    y = state_history(2, k);
    phi = state_history(3, k);
    q1 = state_history(4, k);
    q2 = state_history(5, k);

    % Rotate and translate spacecraft
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
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


% Plot End-Effector Error Over Time
figure;
plot(time(1:size(x_tilde_history, 2)), vecnorm(x_tilde_history, 2, 1), 'k', 'LineWidth', 1.5);
title('End-Effector Position Error Over Time');
xlabel('Time (s)');
ylabel('Error (m)');
grid on;

% Plot Torque Inputs Over Time
figure;
plot(time(1:size(tau_history, 2)), tau_history(1, :), 'r', 'LineWidth', 1.5);
hold on;
plot(time(1:size(tau_history, 2)), tau_history(2, :), 'b', 'LineWidth', 1.5);
title('Torque Inputs Over Time');
xlabel('Time (s)');
ylabel('Torque (N·m)');
legend('Joint 1 Torque', 'Joint 2 Torque');
grid on;


%
%
%
%


function J = mpc_cost_coupled(u, state, horizon, target, Q, R, thruster_penalty, l1, l2, dt)
    % Initialize cost
    J = 0;
    
    % State prediction over the horizon
    state_horizon = state;

    for t = 1:horizon
        % Extract joint torques for this time step
        tau = u(2 * (t - 1) + 1:2 * t);

        % Extract current states
        q1 = state_horizon(4);
        q2 = state_horizon(5);
        dq1 = state_horizon(9);
        dq2 = state_horizon(10);

        % Forward kinematics for the end-effector
        x_end = l1 * cos(q1) + l2 * cos(q1 + q2);
        y_end = l1 * sin(q1) + l2 * sin(q1 + q2);
        error = [x_end; y_end] - target;

        % Increment cost (state error + input cost + thruster penalty)
        J = J + error' * Q * error + tau' * R * tau;

        % Dynamics of the manipulator
        H11 = 5 * 0.5^2 + 5 * (l1^2 + l2^2) + 2 * 5 * l1 * l2 * cos(q2);
        H22 = 5 * l2^2;
        H12 = 5 * l1 * l2 * cos(q2);
        H = [H11, H12; H12, H22];
        h = -5 * l1 * l2 * sin(q2);
        C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];
        ddq = H \ (tau - C * [dq1; dq2]);

        % Reaction torque affecting spacecraft
        tau_reaction = sum(tau);
        domega_dt = tau_reaction / 10.0; % Inertia of the satellite

        % Update joint velocities and positions
        dq1 = dq1 + ddq(1) * dt;
        dq2 = dq2 + ddq(2) * dt;
        q1 = q1 + dq1 * dt;
        q2 = q2 + dq2 * dt;

        % Update spacecraft angular velocity and orientation
        omega = state_horizon(8) + domega_dt * dt;
        phi = state_horizon(3) + omega * dt;

        % Thruster cost for spacecraft position adjustment
        position_error = target - [x_end; y_end];
        thruster_cost = thruster_penalty * norm(position_error);

        % Update cost with thruster penalty
        J = J + thruster_cost;

        % Update the state vector for the next step
        state_horizon(3) = phi; % Orientation
        state_horizon(4) = q1; % Joint angle q1
        state_horizon(5) = q2; % Joint angle q2
        state_horizon(8) = omega; % Angular velocity
        state_horizon(9) = dq1; % Joint velocity q1
        state_horizon(10) = dq2; % Joint velocity q2
    end
end


