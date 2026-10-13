%% ATTITUDE_INITIALS
% Spacecraft attitude subsystem initialization
%
% Defines separate LEADER and FOLLOWER attitude parameters.
% All quantities use SI units.
%
% Quaternion convention used here:
% scalar-first --> q = [q0; q1; q2; q3]
%
% Attitude reference:
% LVLH
%
% Body axes:
% x, y, z with z along the long dimension of the 3U CubeSat.


%% ========================================================================
%  COMMON SPACECRAFT GEOMETRY
%  ========================================================================

body_dimensions = [0.10; 0.10; 0.30];    % [m]

n_panels = 2;

panel_dimensions = [0.30; 0.30];         % [m]
panel_area = prod(panel_dimensions);      % [m^2]

attitude_reference   = 'LVLH';
quaternion_convention = 'scalar-first';


%% ========================================================================
%  LEADER SPACECRAFT
%  ========================================================================

%% Mass
m_sat_leader = 2;                       % [kg]

%% Initial attitude
% Quaternion:
% q = [q0; q1; q2; q3]
%
% [1;0;0;0] means zero attitude error relative to LVLH.

q0_leader = [1;
             0;
             0;
             0];

%% Initial angular velocity
omega0_leader = [0;
                 0;
                 0];                   %PQR    [rad/s]

%% Inertia tensor
% About spacecraft center of mass and expressed in body axes.

J_leader = [0.0078385   0          0;
            0          0.0057552   0;
            0          0          0.0095833];  % [kg*m^2]


%% ========================================================================
%  FOLLOWER SPACECRAFT
%  ========================================================================

%% Mass
m_sat_follower = 2;                     % [kg]

%% Initial attitude
%
% Example:
% q0_follower = [1;0;0;0]
%
% gives the same initial LVLH alignment as the leader.

q0_follower = [1;
               0;
               0;
               0];

%% Initial angular velocity
omega0_follower = [0;
                   0;
                   0];                    % [rad/s]

%% Inertia tensor

J_follower = [0.0078385   0          0;
              0           0.0057552  0;
              0           0          0.0095833]; % [kg*m^2]


%% ========================================================================
%  REACTION-WHEEL TORQUE LIMITS
%  ========================================================================

% Independent symmetric torque saturation:
%
% -Mmax <= M <= +Mmax
%
% for each body axis.

Mmax_leader = [3.2e-3;
               3.2e-3;
               3.2e-3];                  % [N*m]

Mmax_follower = [3.2e-3;
                 3.2e-3;
                 3.2e-3];                % [N*m]


%% ========================================================================
%  LEADER PD CONTROLLER
%  ========================================================================

% Desired closed-loop characteristics

zeta_leader = 0.8;                        % damping ratio
wn_leader   = 0.15;                       % natural frequency [rad/s]

% Principal moments of inertia

Jx_leader = J_leader(1,1);
Jy_leader = J_leader(2,2);
Jz_leader = J_leader(3,3);

% For each axis:
%
% J*theta_ddot + Kd*theta_dot + Kp*theta = 0
%
% Comparing with
%
% theta_ddot + 2*zeta*wn*theta_dot + wn^2*theta = 0
%
% gives
%
% Kp = J*wn^2
% Kd = 2*zeta*J*wn

Kp_leader = [Jx_leader;
             Jy_leader;
             Jz_leader] * wn_leader^2;

Kd_leader = 2*zeta_leader*wn_leader * ...
            [Jx_leader;
             Jy_leader;
             Jz_leader];


%% ========================================================================
%  FOLLOWER PD CONTROLLER
%  ========================================================================

% These can be different from the leader if desired.

zeta_follower = 0.8;                      % damping ratio
wn_follower   = 0.15;                     % natural frequency [rad/s]

Jx_follower = J_follower(1,1);
Jy_follower = J_follower(2,2);
Jz_follower = J_follower(3,3);

Kp_follower = [Jx_follower;
               Jy_follower;
               Jz_follower] * wn_follower^2;

Kd_follower = 2*zeta_follower*wn_follower * ...
              [Jx_follower;
               Jy_follower;
               Jz_follower];


%% ========================================================================
%  NORMALIZE INITIAL QUATERNIONS
%  ========================================================================

q0_leader_norm = norm(q0_leader);

assert(q0_leader_norm > 0, ...
    'Leader initial quaternion must have nonzero norm.');

if abs(q0_leader_norm - 1) > 1e-12
    warning('Leader quaternion was normalized.');
end

q0_leader = q0_leader / q0_leader_norm;


q0_follower_norm = norm(q0_follower);

assert(q0_follower_norm > 0, ...
    'Follower initial quaternion must have nonzero norm.');

if abs(q0_follower_norm - 1) > 1e-12
    warning('Follower quaternion was normalized.');
end

q0_follower = q0_follower / q0_follower_norm;


%% ========================================================================
%  PARAMETER VALIDATION
%  ========================================================================

% ----- Leader -----

assert(isequal(size(q0_leader),[4 1]), ...
    'q0_leader must be 4x1.');

assert(isequal(size(omega0_leader),[3 1]), ...
    'omega0_leader must be 3x1.');

assert(isequal(size(J_leader),[3 3]), ...
    'J_leader must be 3x3.');

assert(norm(J_leader-J_leader.','fro') < 1e-12, ...
    'J_leader must be symmetric.');

[~,chol_flag_leader] = chol(J_leader);

assert(chol_flag_leader == 0, ...
    'J_leader must be positive definite.');


% ----- Follower -----

assert(isequal(size(q0_follower),[4 1]), ...
    'q0_follower must be 4x1.');

assert(isequal(size(omega0_follower),[3 1]), ...
    'omega0_follower must be 3x1.');

assert(isequal(size(J_follower),[3 3]), ...
    'J_follower must be 3x3.');

assert(norm(J_follower-J_follower.','fro') < 1e-12, ...
    'J_follower must be symmetric.');

[~,chol_flag_follower] = chol(J_follower);

assert(chol_flag_follower == 0, ...
    'J_follower must be positive definite.');


%% ========================================================================
%  DISPLAY CONTROLLER PARAMETERS
%  ========================================================================

fprintf('\n====================================================\n');
fprintf('            ATTITUDE INITIALIZATION\n');
fprintf('====================================================\n');


fprintf('\n--- LEADER ---\n');

fprintf('Initial quaternion:\n');
disp(q0_leader);

fprintf('Initial angular velocity [rad/s]:\n');
disp(omega0_leader);

fprintf('Kp roll  = %.8f N*m/rad\n', Kp_leader(1));
fprintf('Kp pitch = %.8f N*m/rad\n', Kp_leader(2));
fprintf('Kp yaw   = %.8f N*m/rad\n', Kp_leader(3));

fprintf('Kd roll  = %.8f N*m*s/rad\n', Kd_leader(1));
fprintf('Kd pitch = %.8f N*m*s/rad\n', Kd_leader(2));
fprintf('Kd yaw   = %.8f N*m*s/rad\n', Kd_leader(3));


fprintf('\n--- FOLLOWER ---\n');

fprintf('Initial quaternion:\n');
disp(q0_follower);

fprintf('Initial angular velocity [rad/s]:\n');
disp(omega0_follower);

fprintf('Kp roll  = %.8f N*m/rad\n', Kp_follower(1));
fprintf('Kp pitch = %.8f N*m/rad\n', Kp_follower(2));
fprintf('Kp yaw   = %.8f N*m/rad\n', Kp_follower(3));

fprintf('Kd roll  = %.8f N*m*s/rad\n', Kd_follower(1));
fprintf('Kd pitch = %.8f N*m*s/rad\n', Kd_follower(2));
fprintf('Kd yaw   = %.8f N*m*s/rad\n', Kd_follower(3));


%% Remove temporary variables

clear q0_leader_norm q0_follower_norm
clear chol_flag_leader chol_flag_follower