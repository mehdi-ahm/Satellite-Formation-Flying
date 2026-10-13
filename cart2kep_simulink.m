function [a,e,i,omega,Omega,nu,M] = cart2kep_simulink(r_vec,v_vec,mu)
% CART2KEP_SIMULINK  Simulink/code-generation friendly Cartesian -> Keplerian.
%
% Inputs
%   r_vec : 3-element inertial position vector
%   v_vec : 3-element inertial velocity vector
%   mu    : gravitational parameter
%
% IMPORTANT:
%   Units must be consistent. For the current thesis model use:
%       r_vec [m], v_vec [m/s], mu [m^3/s^2]
%
% Outputs
%   a      semi-major axis
%   e      eccentricity
%   i      inclination [rad]
%   omega  argument of perigee [rad]
%   Omega  RAAN [rad]
%   nu     true anomaly [rad]
%   M      mean anomaly [rad]
%
% Output order is intentionally:
%   [a,e,i,omega,Omega,nu,M]
% so the first six values can be sent directly to osc2mean:
%   [a;e;i;omega;Omega;nu]
%
%#codegen

% -------------------------------------------------------------------------
% Force fixed 1x3 internal vectors WITHOUT r_vec(:)' or variable-size reshape.
% This is friendly to Simulink dimension propagation whether the input port is
% configured as a 3x1 or 1x3 three-element vector.
% -------------------------------------------------------------------------
r3 = [r_vec(1), r_vec(2), r_vec(3)];
v3 = [v_vec(1), v_vec(2), v_vec(3)];

% Magnitudes and radial velocity
rmag = sqrt(r3(1)^2 + r3(2)^2 + r3(3)^2);
v2   = v3(1)^2 + v3(2)^2 + v3(3)^2;

% Specific angular momentum h = r x v
hvec = [ ...
    r3(2)*v3(3) - r3(3)*v3(2), ...
    r3(3)*v3(1) - r3(1)*v3(3), ...
    r3(1)*v3(2) - r3(2)*v3(1)];

hmag = sqrt(hvec(1)^2 + hvec(2)^2 + hvec(3)^2);

% Inclination: atan2 form is numerically safer than acos(hz/h)
i = atan2(sqrt(hvec(1)^2 + hvec(2)^2), hvec(3));

% Node vector n = k x h = [-hy, hx, 0]
nvec = [-hvec(2), hvec(1), 0.0];
nmag = sqrt(nvec(1)^2 + nvec(2)^2);

% Eccentricity vector e = (v x h)/mu - r/|r|
vxh = [ ...
    v3(2)*hvec(3) - v3(3)*hvec(2), ...
    v3(3)*hvec(1) - v3(1)*hvec(3), ...
    v3(1)*hvec(2) - v3(2)*hvec(1)];

evec = vxh/mu - r3/rmag;
e = sqrt(evec(1)^2 + evec(2)^2 + evec(3)^2);

% Semi-major axis from vis-viva
a = 1.0/(2.0/rmag - v2/mu);

% Small tolerance only for singular geometry handling.
% The current thesis baseline (e=0.001, i=20 deg) is far from these limits.
tol = 1.0e-12;

% -------------------------------------------------------------------------
% RAAN
% -------------------------------------------------------------------------
if nmag > tol
    Omega = wrap2pi_local(atan2(nvec(2),nvec(1)));
else
    % Equatorial orbit: RAAN is undefined. Use zero as a nonsingular
    % convention so downstream code receives a finite value.
    Omega = 0.0;
end

% -------------------------------------------------------------------------
% Argument of perigee
% Use atan2(sin,cos) rather than acos + quadrant branching.
% -------------------------------------------------------------------------
if (e > tol) && (nmag > tol)
    c_omega = (nvec(1)*evec(1) + nvec(2)*evec(2) + nvec(3)*evec(3)) ...
              /(nmag*e);

    nxev = [ ...
        nvec(2)*evec(3) - nvec(3)*evec(2), ...
        nvec(3)*evec(1) - nvec(1)*evec(3), ...
        nvec(1)*evec(2) - nvec(2)*evec(1)];

    s_omega = (nxev(1)*hvec(1) + nxev(2)*hvec(2) + nxev(3)*hvec(3)) ...
              /(nmag*e*hmag);

    omega = wrap2pi_local(atan2(s_omega,c_omega));

elseif e > tol
    % Equatorial eccentric orbit: longitude of periapsis is used with
    % Omega = 0, therefore omega carries that angle.
    omega = wrap2pi_local(atan2(evec(2),evec(1)));

else
    % Circular orbit: periapsis direction is undefined.
    omega = 0.0;
end

% -------------------------------------------------------------------------
% True anomaly
% -------------------------------------------------------------------------
if e > tol
    c_nu = (evec(1)*r3(1) + evec(2)*r3(2) + evec(3)*r3(3))/(e*rmag);

    exr = [ ...
        evec(2)*r3(3) - evec(3)*r3(2), ...
        evec(3)*r3(1) - evec(1)*r3(3), ...
        evec(1)*r3(2) - evec(2)*r3(1)];

    s_nu = (exr(1)*hvec(1) + exr(2)*hvec(2) + exr(3)*hvec(3)) ...
           /(e*rmag*hmag);

    nu = wrap2pi_local(atan2(s_nu,c_nu));

elseif nmag > tol
    % Circular inclined orbit: use argument of latitude as nu while omega=0.
    c_u = (nvec(1)*r3(1) + nvec(2)*r3(2) + nvec(3)*r3(3))/(nmag*rmag);

    nxr = [ ...
        nvec(2)*r3(3) - nvec(3)*r3(2), ...
        nvec(3)*r3(1) - nvec(1)*r3(3), ...
        nvec(1)*r3(2) - nvec(2)*r3(1)];

    s_u = (nxr(1)*hvec(1) + nxr(2)*hvec(2) + nxr(3)*hvec(3)) ...
          /(nmag*rmag*hmag);

    nu = wrap2pi_local(atan2(s_u,c_u));

else
    % Circular equatorial orbit: use true longitude as nu.
    nu = wrap2pi_local(atan2(r3(2),r3(1)));
end

% -------------------------------------------------------------------------
% Mean anomaly for elliptic orbit
% For e=0 this correctly reduces to M = nu.
% -------------------------------------------------------------------------
root_term = max(0.0,1.0-e^2);
E = atan2(sqrt(root_term)*sin(nu), e + cos(nu));
E = wrap2pi_local(E);
M = wrap2pi_local(E - e*sin(E));

end


function ang = wrap2pi_local(ang)
% Fixed-scalar wrap into [0,2*pi)
ang = mod(ang,2.0*pi);
if ang < 0.0
    ang = ang + 2.0*pi;
end
end
