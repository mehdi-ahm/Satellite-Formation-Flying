function [Mb, angle_error, Mb_raw] = attitude_pd_controller( ...
    theta_cmd, attitude_zyx, pqr, Kp, Kd, Mmax)
%ATTITUDE_PD_CONTROLLER Spacecraft attitude PD control with torque limits.
%#codegen
%
% Suitable for a Simulink MATLAB Function block. Configure theta_cmd as
% a scalar and all other inputs as fixed-size [3 1] real double vectors.
% All outputs are fixed-size [3 1] column vectors.
%
% INPUTS
%   theta_cmd    Desired pitch [rad]; desired roll and yaw are zero.
%   attitude_zyx Measured ZYX Euler angles [yaw; pitch; roll] [rad].
%   pqr          Body angular rates [p; q; r] [rad/s].
%   Kp           Proportional gains [roll; pitch; yaw] [N*m/rad].
%   Kd           Derivative gains [roll; pitch; yaw] [N*m*s/rad].
%   Mmax         Nonnegative torque limits [roll; pitch; yaw] [N*m].
%
% OUTPUTS (body-axis order: roll/x, pitch/y, yaw/z)
%   Mb           Saturated commanded body torque [N*m].
%   angle_error  Wrapped measured-minus-command Euler errors [rad].
%   Mb_raw       Unsaturated commanded body torque [N*m].
%
% ASSUMPTIONS
%   Inputs are finite; gains and limits are nonnegative. Body rates are
%   measured relative to the attitude reference frame (e.g. LVLH), and
%   desired relative body rates are zero. No command-rate feedforward is
%   included. Body rates are not generally equal to Euler angle rates.
%   This componentwise Euler PD law is intended for operation away from
%   the ZYX Euler singularity at pitch = +/-pi/2; it is not a global
%   attitude controller. Wrapping each error does not remove that
%   singularity.

% Allocate fixed-size column outputs for MATLAB Coder / Simulink.
Mb = zeros(3, 1);
angle_error = zeros(3, 1);
Mb_raw = zeros(3, 1);

% Reorder [yaw; pitch; roll] into body-axis [roll; pitch; yaw],
% then subtract the command [0; theta_cmd; 0]. Explicit indexing
% avoids implicit expansion and preserves fixed-size column handling.
error_unwrapped = zeros(3, 1);
error_unwrapped(1) = attitude_zyx(3);
error_unwrapped(2) = attitude_zyx(2) - theta_cmd;
error_unwrapped(3) = attitude_zyx(1);

for axis = 1:3
    % Wrap measured-minus-command error to [-pi, pi].
    angle_error(axis) = atan2(sin(error_unwrapped(axis)), ...
                             cos(error_unwrapped(axis)));

    % Componentwise equivalent of -Kp.*angle_error - Kd.*pqr.
    Mb_raw(axis) = -Kp(axis) * angle_error(axis) ...
                  -Kd(axis) * pqr(axis);

    % Independently saturate each body torque to +/-Mmax(axis).
    Mb(axis) = min(max(Mb_raw(axis), -Mmax(axis)), Mmax(axis));
end
end
