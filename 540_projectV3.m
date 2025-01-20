clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
m1 = 5; m2 = 5;           % Link masses
Lambda = diag([10, 10]);  % Gain matrix for sliding variable
Kp_end_effector = 100;    % Proportional gain
Kd_end_effector = 50;     % Damping gain
Gamma = 5 * eye(1);       % Adaptation gain matrix
boundary_limit = 10;      % Bounds for the simulation space

% Simulation parameters
t_final = 20; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, phi, q1, q2, vx, vy, omega, dq1, dq2, theta_hat]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0; 1; 1]; % Initial state

% Debug: Verify initial state size
disp(['Initial state size: ', num2str(size(state, 1))]);

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

% Data storage
state_history = zeros(length(state), num_steps); % Preallocate with size of state
x_tilde_history = zeros(2, num_steps); % End-effector error
theta_hat_history = zeros(2, num_steps); % Adaptive parameter estimates

% Debug: Verify preallocated state_history size
disp(['State history preallocation size: ', num2str(size(state_history, 1)), ' rows, ', num2str(size(state_history, 2)), ' columns']);

for i = 1:num_steps
    % Extract current state
    position = state(1:2);     % Satellite position (2x1)
    phi = state(3);            % Satellite orientation (1x1)
    q1 = state(4); q2 = state(5); % Joint angles (1x1 each)
    velocity = state(6:7);     % Satellite velocity (2x1)
    omega = state(8);          % Satellite angular velocity (1x1)
    dq1 = state(9); dq2 = state(10); % Joint angular velocities (1x1 each)
    theta_hat = state(11:12);  % Adaptive parameter estimates (2x1)

    % Debug: Verify sizes of each component
    disp('--- Debugging state components ---');
    disp(['Position size: ', num2str(size(position, 1)), 'x', num2str(size(position, 2))]);
    disp(['Phi size: ', num2str(size(phi, 1)), 'x', num2str(size(phi, 2))]);
    disp(['q1 size: ', num2str(size(q1, 1)), 'x', num2str(size(q1, 2))]);
    disp(['q2 size: ', num2str(size(q2, 1)), 'x', num2str(size(q2, 2))]);
    disp(['Velocity size: ', num2str(size(velocity, 1)), 'x', num2str(size(velocity, 2))]);
    disp(['Omega size: ', num2str(size(omega, 1)), 'x', num2str(size(omega, 2))]);
    disp(['dq1 size: ', num2str(size(dq1, 1)), 'x', num2str(size(dq1, 2))]);
    disp(['dq2 size: ', num2str(size(dq2, 1)), 'x', num2str(size(dq2, 2))]);
    disp(['Theta_hat size: ', num2str(size(theta_hat, 1)), 'x', num2str(size(theta_hat, 2))]);
    
    % Forward kinematics
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    end_effector = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];

    % Check if the end-effector has reached the target
    if norm(end_effector - target_position) < goal_tolerance
        disp('End effector reached the target!');
        state_history = state_history(:, 1:i); % Trim unused steps
        x_tilde_history = x_tilde_history(:, 1:i); % Trim error history
        theta_hat_history = theta_hat_history(:, 1:i); % Trim parameter estimates
        break;
    end

    % Compute position and velocity errors
    x_tilde = end_effector - target_position; % Position error
    dx_tilde = velocity; % Velocity error (simplified here)

    % Sliding variable
    s = dx_tilde + Lambda * x_tilde;

    % Jacobian matrix for the 2-link manipulator
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];

    % Adaptive control law
    tau_manipulator = J' * (-Kp_end_effector * x_tilde - Kd_end_effector * dx_tilde - theta_hat' * s);

     

    % Update parameter estimates (adaptation law)
    dtheta_hat = -Gamma * (s' * dx_tilde);
    theta_hat = theta_hat + dtheta_hat * dt;

    % Manipulator dynamics
    H11 = m1*0.5^2 + m2*(l1^2 + l2^2) + 2*m2*l1*l2*cos(q2) + inertia_satellite;
    H22 = m2*l2^2 + inertia_satellite;
    H12 = m2*l1*l2*cos(q2);
    H = [H11, H12; H12, H22];

    % Coriolis terms
    h = -m2 * l1 * l2 * sin(q2);
    C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];

    % Joint accelerations
    ddq = H \ (tau_manipulator - C * [dq1; dq2]);

    % Update manipulator states
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Update spacecraft dynamics using manipulator reaction torque
    omega = omega + sum(tau_manipulator) / inertia_satellite * dt;
    phi = phi + omega * dt;

    % Ensure theta_hat is a column vector
    theta_hat = theta_hat(:); % Force theta_hat to be a column vector

    % Debug: Verify updated theta_hat size
    disp(['Updated Theta_hat size: ', num2str(size(theta_hat, 1))]);

    % Update state
    state = [position; phi; q1; q2; velocity; omega; dq1; dq2; theta_hat];

    % Debug: Verify state size
    disp(['Updated state size: ', num2str(size(state, 1))]);

    % Check state_history size before assignment
    disp(['State history size: ', num2str(size(state_history, 1))]);
    if size(state, 1) ~= size(state_history, 1)
        error('State vector size (%d) does not match state_history size (%d)', size(state, 1), size(state_history, 1));
    end

    % Assign state to history
    state_history(:, i) = state;  % Store current state

    % Log data
    x_tilde_history(:, i) = x_tilde;
    theta_hat_history(:, i) = theta_hat;
end

% Visualization
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft and 2-Link Manipulator with Reaction Control');
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
h_link2 = plot([0, 0], [0, 0], 'k', 'LineWidth', 2); % Link 2

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


%%

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;       % Mass of the satellite
inertia_satellite = 15.0;    % Moment of inertia of the satellite
l1 = 1; l2 = 1;              % Link lengths
m1 = 5; m2 = 5;              % Link masses
Lambda = diag([5, 5]);       % Sliding gain matrix
Kp_thrust = 10;              % Thruster proportional gain
Kp_orientation = 10;         % Orientation proportional gain
Gamma = diag([5, 5]);        % Adaptation gain matrix
boundary_limit = 10;         % Simulation space bounds

% Simulation Parameters
dt = 0.1;                    % Time step
t_final = 30;                % Simulation duration
num_steps = t_final / dt;    % Number of time steps
time = 0:dt:t_final;         % Time vector

% Initial Conditions
% State: [x, y, phi, q1, q2, vx, vy, omega, dq1, dq2, theta_hat, phi_hat]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0; 1; 1];

% Target position
target_position = [5; 5];    % Target in 2D
goal_tolerance = 0.1;        % Tolerance for reaching the target

% Preallocate for Storage
state_history = zeros(length(state), num_steps);
theta_hat_history = zeros(2, num_steps);

% Main Loop
for step = 1:num_steps
    % Extract current state
    position = state(1:2);         % [x, y] (2x1 vector)
    phi = state(3);                % Orientation (scalar)
    q1 = state(4); q2 = state(5);  % Joint angles (scalars)
    velocity = state(6:7);         % Translational velocity [vx, vy] (2x1 vector)
    omega = state(8);              % Rotational velocity (scalar)
    dq1 = state(9); dq2 = state(10); % Joint angular velocities (scalars)
    theta_hat = state(11);         % Adaptive parameter for joints (scalar)
    phi_hat = state(12);           % Adaptive parameter for spacecraft orientation (scalar)

    % --- Debugging Dimensions ---
    fprintf('Step %d\n', step);
    fprintf('position: %s\n', mat2str(size(position)));
    fprintf('phi: %s\n', mat2str(size(phi)));
    fprintf('q1: %s\n', mat2str(size(q1)));
    fprintf('q2: %s\n', mat2str(size(q2)));
    fprintf('velocity: %s\n', mat2str(size(velocity)));
    fprintf('omega: %s\n', mat2str(size(omega)));
    fprintf('dq1: %s\n', mat2str(size(dq1)));
    fprintf('dq2: %s\n', mat2str(size(dq2)));
    fprintf('theta_hat: %s\n', mat2str(size(theta_hat)));
    fprintf('phi_hat: %s\n', mat2str(size(phi_hat)));

    % Enforce Consistent Dimensions
    position = position(:);    % Ensure column vector
    velocity = velocity(:);    % Ensure column vector
    theta_hat = theta_hat(:);  % Ensure scalar
    phi_hat = phi_hat(:);      % Ensure scalar

    % Forward Kinematics
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    end_effector = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];

    % Check Target Reached
    if norm(end_effector - target_position) < goal_tolerance
        disp('Target reached!');
        state_history = state_history(:, 1:step);
        theta_hat_history = theta_hat_history(:, 1:step);
        break;
    end

    % Errors
    x_tilde = end_effector - target_position;        % Position error
    dx_tilde = velocity;                             % Velocity error

    % Sliding Variable
    s = dx_tilde + Lambda * x_tilde;

    % Jacobian of the Manipulator
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];

    % Adaptive Controller
    % Control torque for the manipulator
    tau_manipulator = -J' * (Kp_orientation * x_tilde + theta_hat * s);

    % Update adaptive parameters (adaptation law)
    scalar_product_joint = s' * dx_tilde;            % Scalar product for joint dynamics
    scalar_product_orientation = s' * omega;        % Scalar product for orientation dynamics

    % Update rules for adaptive parameters
    dtheta_hat = -Gamma(1, 1) * scalar_product_joint; % Joint dynamics
    dphi_hat = -Gamma(2, 2) * scalar_product_orientation; % Orientation dynamics

    % Update adaptive parameters
    theta_hat = theta_hat + dtheta_hat * dt;
    phi_hat = phi_hat + dphi_hat * dt;

    % Dynamics of the Manipulator
    % Inertia matrix
    H11 = m1 * (l1 / 2)^2 + m2 * (l1^2 + (l2 / 2)^2 + 2 * l1 * (l2 / 2) * cos(q2));
    H22 = m2 * (l2 / 2)^2;
    H12 = m2 * l1 * (l2 / 2) * cos(q2);
    H = [H11, H12; H12, H22];

    % Coriolis terms
    h = -m2 * l1 * (l2 / 2) * sin(q2);
    C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];

    % Joint accelerations
    ddq = H \ (tau_manipulator - C * [dq1; dq2]);

    % Update Joint States
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Thruster Control
    thrust_force = -Kp_thrust * x_tilde;             % Translational control
    dposition_dt = velocity;                         % Velocity affects position
    dvelocity_dt = thrust_force / mass_satellite;    % Force affects velocity

    % Update Satellite States
    omega = omega + phi_hat * dt; % Orientation torque from adaptive controller
    phi = phi + omega * dt;

    velocity = velocity + dvelocity_dt * dt;
    position = position + dposition_dt * dt;

    % Ensure scalars are column vectors
    phi = reshape(phi, 1, 1);
    q1 = reshape(q1, 1, 1);
    q2 = reshape(q2, 1, 1);
    omega = reshape(omega, 1, 1);
    dq1 = reshape(dq1, 1, 1);
    dq2 = reshape(dq2, 1, 1);
    theta_hat = reshape(theta_hat, 1, 1);
    phi_hat = reshape(phi_hat, 1, 1);
    
    % Update State
    state = [position; phi; q1; q2; velocity; omega; dq1; dq2; theta_hat; phi_hat];

    % Log Data
    state_history(:, step) = state;
    theta_hat_history(:, step) = [theta_hat; phi_hat];
end

% Visualization of Adaptive Parameters
figure;
plot(time(1:size(theta_hat_history, 2)), theta_hat_history(1, :), 'r', ...
     time(1:size(theta_hat_history, 2)), theta_hat_history(2, :), 'b');
title('Adaptive Parameter Estimates');
xlabel('Time (s)');
ylabel('Theta Hat');
legend('\theta_{hat} for Joints', '\phi_{hat} for Orientation');

%%
clc;
clear;
close all;

% Constants
mass_satellite = 10.0;       % Mass of the satellite
inertia_satellite = 15.0;    % Moment of inertia of the satellite
l1 = 1; l2 = 1;              % Link lengths
m1 = 5; m2 = 5;              % Link masses
Lambda = diag([5, 5]);       % Sliding gain matrix
Kp_thrust = 10;              % Thruster proportional gain
Kp_orientation = 10;         % Orientation proportional gain
Gamma = diag([5, 5]);        % Adaptation gain matrix
boundary_limit = 10;         % Simulation space bounds

% Simulation Parameters
dt = 0.1;                    % Time step
t_final = 30;                % Simulation duration
num_steps = t_final / dt;    % Number of time steps
time = 0:dt:t_final;         % Time vector

% Initial Conditions
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0; 1; 1]; % State: [x, y, phi, q1, q2, vx, vy, omega, dq1, dq2, theta_hat, phi_hat]

% Target position
target_position = [5; 5];    % Target in 2D
goal_tolerance = 0.1;        % Tolerance for reaching the target

% Preallocate for Storage
state_history = zeros(length(state), num_steps);
theta_hat_history = zeros(2, num_steps);

% Main Loop
for step = 1:num_steps
    % Extract current state
    position = state(1:2);
    phi = state(3);
    q1 = state(4); q2 = state(5);
    velocity = state(6:7);
    omega = state(8);
    dq1 = state(9); dq2 = state(10);
    theta_hat = state(11);
    phi_hat = state(12);

    % Ensure all state variables are properly formatted
    position = position(:);    % Ensure column vector
    velocity = velocity(:);    % Ensure column vector
    % Scalars, no need for reshape
    %phi = reshape(phi, 1, 1);
    %q1 = reshape(q1, 1, 1);
    %q2 = reshape(q2, 1, 1);
    %omega = reshape(omega, 1, 1);
    %dq1 = reshape(dq1, 1, 1);
    %dq2 = reshape(dq2, 1, 1);
    %theta_hat = reshape(theta_hat, 1, 1);
    %phi_hat = reshape(phi_hat, 1, 1);

    % Forward Kinematics
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p1 = position + R * [l1 * cos(q1); l1 * sin(q1)];
    end_effector = p1 + R * [l2 * cos(q1 + q2); l2 * sin(q1 + q2)];

    % Check Target Reached
    if norm(end_effector - target_position) < goal_tolerance
        disp('Target reached!');
        state_history = state_history(:, 1:step);
        theta_hat_history = theta_hat_history(:, 1:step);
        break;
    end

    % Errors
    x_tilde = end_effector - target_position;
    dx_tilde = velocity;

    % Sliding Variable
    s = dx_tilde + Lambda * x_tilde;

    % Jacobian of the Manipulator
    J = [
        -l1 * sin(q1) - l2 * sin(q1 + q2), -l2 * sin(q1 + q2);
        l1 * cos(q1) + l2 * cos(q1 + q2),  l2 * cos(q1 + q2)
    ];

    % Adaptive Controller
    tau_manipulator = -J' * (Kp_orientation * x_tilde + theta_hat * s);

    % Update adaptive parameters (adaptation law)
    scalar_product_joint = s' * dx_tilde;
    scalar_product_orientation = s' * omega;

    % Update rules for adaptive parameters
    dtheta_hat = -Gamma(1, 1) * scalar_product_joint;
    dphi_hat = -Gamma(2, 2) * scalar_product_orientation;

    % Update adaptive parameters
    theta_hat = theta_hat + dtheta_hat * dt;
    phi_hat = phi_hat + dphi_hat * dt;

    % Dynamics of the Manipulator
    H11 = m1 * (l1 / 2)^2 + m2 * (l1^2 + (l2 / 2)^2 + 2 * l1 * (l2 / 2) * cos(q2));
    H22 = m2 * (l2 / 2)^2;
    H12 = m2 * l1 * (l2 / 2) * cos(q2);
    H = [H11, H12; H12, H22];

    % Coriolis terms
    h = -m2 * l1 * (l2 / 2) * sin(q2);
    C = [h * dq2, h * (dq1 + dq2); -h * dq1, 0];

    % Joint accelerations
    ddq = H \ (tau_manipulator - C * [dq1; dq2]);

    % Update Joint States
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Thruster Control
    thrust_force = -Kp_thrust * x_tilde;
    dposition_dt = velocity;
    dvelocity_dt = thrust_force / mass_satellite;

    % Update Satellite States
    omega = omega + phi_hat * dt;
    phi = phi + omega * dt;
    velocity = velocity + dvelocity_dt * dt;
    position = position + dposition_dt * dt;

    % Update State
    % Ensure all components are explicitly column vectors or scalars
    position = position(:);    % Force 2x1 column vector
    velocity = velocity(:);    % Force 2x1 column vector
    phi = phi(:);              % Force scalar as column vector (1x1)
    q1 = q1(:);                % Force scalar as column vector (1x1)
    q2 = q2(:);                % Force scalar as column vector (1x1)
    omega = omega(:);          % Force scalar as column vector (1x1)
    dq1 = dq1(:);              % Force scalar as column vector (1x1)
    dq2 = dq2(:);              % Force scalar as column vector (1x1)
    theta_hat = theta_hat(:);  % Ensure column vector (1x1)
    phi_hat = phi_hat(:);      % Ensure column vector (1x1)
    
    % Concatenate all components into the state vector
    state = [position; phi; q1; q2; velocity; omega; dq1; dq2; theta_hat; phi_hat];

    % % Separate history storage and reload state
    % state = [position; phi; q1; q2; velocity; omega; dq1; dq2; theta_hat; phi_hat];

    % Log Data
    state_history(:, step) = state;
    theta_hat_history(:, step) = [theta_hat; phi_hat];
end

% Visualization of Adaptive Parameters
figure;
plot(time(1:size(theta_hat_history, 2)), theta_hat_history(1, :), 'r', ...
     time(1:size(theta_hat_history, 2)), theta_hat_history(2, :), 'b');
title('Adaptive Parameter Estimates');
xlabel('Time (s)');
ylabel('Theta Hat');
legend('\theta_{hat} for Joints', '\phi_{hat} for Orientation');

%% WORKING ADAPTIVE CONTROL !!!!
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
thruster_force = 0.2;     % Thruster boost
boundary_limit = 10;      % Bounds for the simulation space

% Limits
q_limit = deg2rad(150);    % Maximum joint angle limit
dq_limit = deg2rad(50);    % Maximum joint velocity limit
omega_limit = 1;           % Maximum satellite angular velocity

% Simulation parameters
t_final = 100; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Adaptive parameters (initial estimates)
m1_hat = 1.0; % Estimated mass of link 1
m2_hat = 1.0; % Estimated mass of link 2

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

% Data storage for analysis and visualization
state_history = zeros(length(state), num_steps);
x_tilde_history = zeros(2, num_steps); % End-effector error
m1_hat_history = zeros(1, num_steps); % Adaptive parameter m1_hat
m2_hat_history = zeros(1, num_steps); % Adaptive parameter m2_hat

for i = 1:num_steps
    % Store state
    state_history(:, i) = state;
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
end

% Visualization
figure;
hold on;
axis equal;
xlim([-boundary_limit, boundary_limit]);
ylim([-boundary_limit, boundary_limit]);
title('Spacecraft and 2-Link Manipulator with Reaction Control');
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


% Parameter adaptation and error plots
figure;

% % Adaptive parameter estimates (m1_hat and m2_hat)
% subplot(3, 1, 1);
% plot(time(1:length(m1_hat_history)), m1_hat_history, 'r', 'LineWidth', 1.5);
% hold on;
% plot(time(1:length(m2_hat_history)), m2_hat_history, 'b', 'LineWidth', 1.5);
% title('Adaptive Parameters Over Time (Mass Estimates)');
% xlabel('Time (s)');
% ylabel('Mass Estimates');
% legend('m1\_hat', 'm2\_hat');
% grid on;
% 
% % Adaptive parameter estimates for q1 and q2
% subplot(3, 1, 2);
% q1_history = rad2deg(state_history(4, 1:size(m1_hat_history, 2))); % q1 in degrees
% q2_history = rad2deg(state_history(5, 1:size(m2_hat_history, 2))); % q2 in degrees
% plot(time(1:length(q1_history)), q1_history, 'm', 'LineWidth', 1.5);
% hold on;
% plot(time(1:length(q2_history)), q2_history, 'c', 'LineWidth', 1.5);
% title('Adaptive Parameters for Joint Angles Over Time');
% xlabel('Time (s)');
% ylabel('Joint Angles (Degrees)');
% legend('q1', 'q2');
% grid on;
% 
% % End-effector error over time
% subplot(3, 1, 3);
plot(time(1:size(x_tilde_history, 2)), vecnorm(x_tilde_history, 2, 1), 'k', 'LineWidth', 1.5);
title('End-Effector Error Over Time');
xlabel('Time (s)');
ylabel('Error (m)');
grid on;



%% Slightly simplified Adaptive controller (not as robust though)

clc;
clear;
close all;

% Constants
mass_satellite = 10.0;    % Mass of the satellite
inertia_satellite = 10.0; % Moment of inertia of the satellite
l1 = 1; l2 = 1;           % Link lengths
m1_hat = 5; m2_hat = 5;   % Initial adaptive mass estimates
Kp_end_effector = 100;    % End-effector proportional gain
Kd_end_effector = 50;     % End-effector damping gain
gamma = 0.1;              % Adaptation gain
thruster_force = 0.2;     % Thruster boost
boundary_limit = 10;      % Bounds for the simulation space

% Limits
q_limit = deg2rad(150);    % Maximum joint angle limit
dq_limit = deg2rad(50);    % Maximum joint velocity limit
omega_limit = 1;           % Maximum satellite angular velocity

% Simulation parameters
t_final = 200; dt = 0.1;
num_steps = t_final / dt;
time = 0:dt:t_final-dt; % Time vector

% Initial conditions: [x, y, theta, q1, q2, vx, vy, omega, dq1, dq2]
state = [0; 0; 0; deg2rad(30); deg2rad(-45); 0; 0; 0; 0; 0]; % Initial state

% Target position
target_position = [4; 4]; % 2D position
goal_tolerance = 0.1;

% Initialize History Arrays
state_history = zeros(length(state), num_steps);
x_tilde_history = zeros(2, num_steps); % End-effector error history
m1_hat_history = zeros(1, num_steps);  % Adaptive mass estimate for link 1
m2_hat_history = zeros(1, num_steps);  % Adaptive mass estimate for link 2

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
        m1_hat_history = m1_hat_history(1:i); % Trim parameter history
        m2_hat_history = m2_hat_history(1:i); % Trim parameter history
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

    % Adaptive parameter update laws
    m1_hat_dot = -gamma * (x_tilde' * [l1 * cos(q1); 0]);
    m2_hat_dot = -gamma * (x_tilde' * [l2 * cos(q1 + q2); 0]);
    m1_hat = m1_hat + m1_hat_dot * dt;
    m2_hat = m2_hat + m2_hat_dot * dt;

    % Store adaptive parameter values
    m1_hat_history(i) = m1_hat;
    m2_hat_history(i) = m2_hat;

    % PD control for the manipulator in end-effector space
    tau_manipulator = -J' * (Kp_end_effector * x_tilde + Kd_end_effector * dx_tilde);

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

    % Update manipulator state
    dq1 = dq1 + ddq(1) * dt; q1 = q1 + dq1 * dt;
    dq2 = dq2 + ddq(2) * dt; q2 = q2 + dq2 * dt;

    % Update spacecraft dynamics
    domega_dt = sum(tau_manipulator) / inertia_satellite; % Reaction torque
    omega = omega + domega_dt * dt;
    phi = phi + omega * dt;

    % Thruster control
    position_error = target_position - position;
    thrust = thruster_force * position_error / max(norm(position_error), 1e-6);

    % Compute spacecraft translational dynamics
    dposition_dt = velocity;
    dvelocity_dt = thrust / mass_satellite;
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
title('Spacecraft and 2-Link Manipulator with Reaction Control');
xlabel('X Position (m)');
ylabel('Y Position (m)');

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

% Adaptive parameter and error plots
figure;

% Adaptive parameter estimates
subplot(2, 1, 1);
plot(time(1:length(m1_hat_history)), m1_hat_history, 'r', 'LineWidth', 1.5);
hold on;
plot(time(1:length(m2_hat_history)), m2_hat_history, 'b', 'LineWidth', 1.5);
title('Adaptive Parameters Over Time');
xlabel('Time (s)');
ylabel('Parameter Value');
legend('m1\_hat', 'm2\_hat');
grid on;

% End-effector error
subplot(2, 1, 2);
plot(time(1:size(x_tilde_history, 2)), vecnorm(x_tilde_history, 2, 1), 'k', 'LineWidth', 1.5);
title('End-Effector Error Over Time');
xlabel('Time (s)');
ylabel('Error (m)');
grid on;

