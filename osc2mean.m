function mean_elements = osc2mean(oeosc, j2, req)
% OSC2MEAN  Convert osculating classical orbital elements to Brouwer-Lyddane
%           "mean" elements, using first-order J2 perturbation theory.
%
% Same algorithm as Fanebust's original osc2mean.m (Brouwer-Lyddane
% short-period correction, "Orbital Mechanics with MATLAB" / D. Eagle).
% Only variable names have been changed for readability — every formula,
% branch condition, and iteration count is identical to the original.
%
% Two branches:
%   - e < 0.02  : near-circular orbits. Argument of perigee is undefined
%                 as e -> 0, so this branch instead uses the mean
%                 longitude (mean_lon = M + argp) and the eccentricity
%                 vector components (ecc_x = e*cos(argp), ecc_y = e*sin(argp)),
%                 which stay well-defined at e = 0.
%   - e >= 0.02 : general eccentricity. Works directly with the classical
%                 elements and the true anomaly.
%
% Both branches iterate the correction 5 times: each pass recomputes the
% correction terms using the previous pass's "mean" estimate.
%
% INPUT
%  oeosc(1) = osculating semimajor axis           (same length unit as req)
%  oeosc(2) = osculating eccentricity                    [-]   0 <= e < 1
%  oeosc(3) = osculating inclination                    [rad]  0 <= i <= pi
%  oeosc(4) = osculating argument of perigee            [rad]  0 <= 2*pi
%  oeosc(5) = osculating right ascension of asc. node    [rad]  0 <= 2*pi
%  oeosc(6) = osculating true anomaly                    [rad]  0 <= 2*pi
%  j2       = J2 zonal harmonic coefficient               [-]
%  req      = equatorial (or reference) radius of the primary body,
%             same length unit as oeosc(1)
%
% OUTPUT
%  mean_elements(1) = mean semimajor axis     (same length unit as req)
%  mean_elements(2) = mean eccentricity                    [-]
%  mean_elements(3) = mean inclination                    [rad]
%  mean_elements(4) = mean argument of perigee            [rad]  0 <= 2*pi
%  mean_elements(5) = mean RAAN                           [rad]  0 <= 2*pi
%  mean_elements(6) = mean mean-anomaly                   [rad]  0 <= 2*pi
%  mean_elements(7) = mean true anomaly (from mean_elements(6) via Kepler's eqn)

%% ------------------------------------------------------------------
%  Step 0 — unpack inputs, convert the given true anomaly to mean anomaly
%  ------------------------------------------------------------------
two_pi = 2.0 * pi;

osc = zeros(6, 1);
for idx = 1:1:6
    osc(idx) = oeosc(idx);
end

% true anomaly -> eccentric anomaly (via a 4-quadrant arctangent built
% from sin/cos of the eccentric anomaly)
sinE_term = sin(osc(6)) * sqrt(1.0 - osc(2) * osc(2));
cosE_term = osc(2) + cos(osc(6));
ecc_anom  = four_quadrant_atan(sinE_term, cosE_term);   % eccentric anomaly

% eccentric anomaly -> mean anomaly (Kepler's equation, direct form)
osc(6) = mod(ecc_anom - osc(2) * sin(ecc_anom), two_pi);

% Short, readable names for the osculating elements ("_osc" suffix).
a_osc    = osc(1);   % osculating semimajor axis
e_osc    = osc(2);   % osculating eccentricity
i_osc    = osc(3);   % osculating inclination
argp_osc = osc(4);   % osculating argument of perigee
raan_osc = osc(5);   % osculating RAAN
M_osc    = osc(6);   % osculating mean anomaly (just computed above)

% Inclination-dependent constants used throughout the J2 correction
% terms, from the standard J2 potential expansion:
%   incl_A = 1/3 - (1/2)sin^2(i),   incl_B = (1/2)sin^2(i)
incl_A = 1.0 / 3.0 - 0.5 * sin(i_osc)^2;
incl_B = 0.5 * sin(i_osc)^2;


if (e_osc < 0.02)
    %% ------------------------------------------------------------------
    %  BRANCH 1 — near-circular orbit (e < 0.01)
    %  Uses mean longitude and eccentricity-vector components instead of
    %  argument of perigee directly, since the latter is ill-defined at
    %  e = 0.
    %  ------------------------------------------------------------------

    mean_lon_osc = mod(M_osc + argp_osc, two_pi);   % osculating mean longitude
    ecc_x_osc    = e_osc * cos(argp_osc);           % eccentricity vector, x
    ecc_y_osc    = e_osc * sin(argp_osc);           % eccentricity vector, y

    % Trig shorthands of the mean longitude, refreshed each iteration.
    sin_lon  = sin(mean_lon_osc);   cos_lon  = cos(mean_lon_osc);
    sin_2lon = sin(2*mean_lon_osc); cos_2lon = cos(2*mean_lon_osc);
    sin_3lon = sin(3*mean_lon_osc); cos_3lon = cos(3*mean_lon_osc);
    sin_4lon = sin(4*mean_lon_osc); cos_4lon = cos(4*mean_lon_osc); %#ok<NASGU>

    sin_2incl = sin(2 * i_osc);
    cos_incl  = cos(i_osc); %#ok<NASGU>

    % Initialize the "mean" estimates with the osculating values —
    % refined by the fixed-point loop below.
    a_mean        = a_osc;
    i_mean        = i_osc;
    raan_mean     = raan_osc;
    mean_lon_mean = mean_lon_osc;
    ecc_x_mean    = ecc_x_osc;
    ecc_y_mean    = ecc_y_osc;

    % ---- fixed-point iteration: 5 passes of the short-period correction ----
    for iter = 1:1:5

        % --- semimajor axis short-period correction, then subtract ---
        da_sp = 3 * j2 * req * req / a_mean * (incl_B * cos_2lon + (1 - 3.5 * incl_B) ...
            * ecc_x_mean * cos_lon + (1 - 2.5 * incl_B) * ecc_y_mean * sin_lon + 3.5 * incl_B ...
            * (ecc_x_mean * cos_3lon + ecc_y_mean * sin_3lon));
        a_mean = a_osc - da_sp;

        % --- inclination short-period correction, then subtract ---
        di_sp = 3 * j2 / 8 * req * req / a_mean^2 * sin_2incl * (cos_2lon - ecc_x_mean * cos_lon ...
            + ecc_y_mean * sin_lon + 7/3 * (ecc_x_mean * cos_3lon + ecc_y_mean * sin_3lon));
        i_mean = i_osc - di_sp;

        % Refresh inclination-dependent constants with the updated i_mean
        cos_incl  = cos(i_mean); %#ok<NASGU>
        sin_2incl = sin(2 * i_mean);
        incl_B    = 0.5 * sin(i_mean)^2;

        % --- RAAN short-period correction, then subtract ---
        draan_sp = 1.5 * j2 * req * req / a_mean^2 * cos(i_mean) ...
            * (0.5 * sin_2lon - 3.5 * ecc_x_mean * sin_lon + 2.5 * ecc_y_mean * cos_lon ...
            + 7/6 * (ecc_x_mean * sin_3lon - ecc_y_mean * cos_3lon));
        raan_mean = raan_osc - draan_sp;

        % --- mean-longitude short-period correction, then subtract ---
        dlon_sp = 1.5 * j2 * req * req / (a_mean * a_mean) ...
            * (-0.5 * (1 - 5 * incl_B) * sin_2lon + (7 - 77/4 * incl_B) ...
            * ecc_x_mean * sin_lon - (6 - 55/4 * incl_B) * ecc_y_mean * cos_lon - (7/6 - 77/12 * incl_B) ...
            * (ecc_x_mean * sin_3lon - ecc_y_mean * cos_3lon));
        mean_lon_mean = mean_lon_osc - dlon_sp;

        % Refresh trig shorthands using the updated mean longitude, for
        % use in the next two corrections (and the next loop pass).
        sin_lon  = sin(mean_lon_mean);   cos_lon  = cos(mean_lon_mean);
        sin_2lon = sin(2*mean_lon_mean); cos_2lon = cos(2*mean_lon_mean);
        sin_3lon = sin(3*mean_lon_mean); cos_3lon = cos(3*mean_lon_mean);
        sin_4lon = sin(4*mean_lon_mean); cos_4lon = cos(4*mean_lon_mean);

        % --- eccentricity-vector x-component correction, then subtract ---
        dex_sp = 1.5 * j2 * req * req / a_mean^2 * ((1 - 2.5 * incl_B) * cos_lon ...
            + 7/6 * incl_B * cos_3lon + (1.5 - 5 * incl_B) * ecc_x_mean * cos_2lon + (2 - 3 * incl_B) ...
            * ecc_y_mean * sin_2lon + 17/4 * incl_B * (ecc_x_mean * cos_4lon + ecc_y_mean * sin_4lon));
        ecc_x_mean = ecc_x_osc - dex_sp;

        % --- eccentricity-vector y-component correction, then subtract ---
        dey_sp = 1.5 * j2 * req * req / a_mean^2 * ((1 - 3.5 * incl_B) * sin_lon ...
            + 7/6 * incl_B * sin_3lon + (1 - 6 * incl_B) * ecc_x_mean * sin_2lon - (1.5 - 4 * incl_B) ...
            * ecc_y_mean * cos_2lon + 17/4 * incl_B * (ecc_x_mean * sin_4lon - ecc_y_mean * cos_4lon));
        ecc_y_mean = ecc_y_osc - dey_sp;
    end
    % end of the 5-pass fixed-point loop

    % Recover mean eccentricity and mean argument of perigee from the
    % refined eccentricity-vector components.
    e_mean = sqrt(ecc_y_mean^2 + ecc_x_mean^2);

    argp_mean = 0;
    if (e_mean > 1.0e-8)
        argp_mean = atan2(ecc_y_mean, ecc_x_mean);
    end

    % "mean" mean anomaly = mean longitude minus mean argument of perigee
    M_mean = mod(mean_lon_mean - argp_mean, two_pi);

else
    %% ------------------------------------------------------------------
    %  BRANCH 2 — general eccentricity (e >= 0.01)
    %  Works directly with classical elements; u_mean = argument of
    %  latitude (argp + true anomaly), r_mean = orbit radius there.
    %  ------------------------------------------------------------------

    p_mean    = a_osc * (1 - e_osc^2);   % osculating semi-latus rectum
    a_mean    = a_osc;
    e_mean    = e_osc;
    i_mean    = i_osc;
    argp_mean = argp_osc;
    raan_mean = raan_osc;
    M_mean    = M_osc;

    % Mean anomaly -> true anomaly via Kepler's equation
    [~, nu_mean] = solve_kepler_equation(M_mean, e_mean);

    u_mean = mod(argp_mean + nu_mean, two_pi);        % argument of latitude
    r_mean = p_mean / (1 + e_mean * cos(nu_mean));    % instantaneous orbit radius

    % ---- fixed-point iteration: 5 passes of the short-period correction ----
    for iter = 1:1:5

        % --- semimajor axis correction ---
        da_sp = 3*j2*req*req/a_mean*((a_mean/r_mean)^3*(incl_A+incl_B*cos(2*u_mean)) ...
            -incl_A*(1-e_mean^2)^(-1.5));
        a_mean = a_osc - da_sp;

        % --- inclination correction ---
        di_sp = 3/8*j2*req^2/p_mean^2*sin(2*i_mean)*(cos(2*u_mean) ...
            +e_mean*cos(nu_mean+2*argp_mean)+1/3*e_mean*cos(3*nu_mean+2*argp_mean));
        i_mean = i_osc - di_sp;

        % Refresh inclination-dependent constants
        incl_A = 1/3 - 0.5*sin(i_mean)^2;
        incl_B = 0.5*sin(i_mean)^2;

        % --- eccentricity correction ---
        de_sp = 1.5*j2*req*req/a_mean^2*((1-e_mean^2)/e_mean*((a_mean/r_mean)^3 ...
            *(incl_A+incl_B*cos(2*u_mean))-incl_A*(1-e_mean^2)^(-1.5))-incl_B/(e_mean*(1-e_mean^2)) ...
            *(cos(2*u_mean)+e_mean*cos(nu_mean+2*argp_mean)+e_mean*cos(3*nu_mean+2.*argp_mean)/3));
        e_mean = e_osc - de_sp;

        % Update the semi-latus rectum and radius with the refined a, e
        p_mean = a_mean * (1.0 - e_mean^2);
        r_mean = p_mean / (1.0 + e_mean * cos(nu_mean));

        nu_mean = mod(nu_mean, 2.0 * pi);
        M_mean  = mod(M_mean, 2.0 * pi);

        % Equation-of-center (true anomaly minus mean anomaly), with a
        % safeguard: near nu_mean or M_mean = 0 or pi, this difference is
        % numerically noisy, so it's forced to zero there instead.
        if ((abs(nu_mean-pi) <= 1.0e-06) || (abs(M_mean-pi) <= 1.0e-06) || ...
                (abs(nu_mean) <= 1.0e-06) || (abs(M_mean) <= 1.0e-06))
            eq_center = 0;
        else
            eq_center = nu_mean - M_mean;
        end

        % --- RAAN correction ---
        draan_sp = -1.5*j2*(req/p_mean)^2*cos(i_mean)*(eq_center+e_mean*sin(nu_mean) ...
            -0.5*sin(2*u_mean)-0.5*e_mean*sin(nu_mean+2*argp_mean)-1/6*e_mean*sin(3*nu_mean+2*argp_mean));
        raan_mean = raan_osc - draan_sp;

        % --- argument of perigee correction ---
        dargp_sp = 1.5*j2*(req/p_mean)^2*((2-5*incl_B)*(eq_center+e_mean*sin(nu_mean)) ...
            +(1-3*incl_B)*((1-0.25*e_mean*e_mean)*sin(nu_mean)/e_mean+0.5*sin(2*nu_mean) ...
            +e_mean*sin(3*nu_mean)/12)-(0.5*incl_B+(0.5-15/8*incl_B)*e_mean^2)/e_mean*sin(nu_mean+2*argp_mean) ...
            +e_mean/8*incl_B*sin(nu_mean-2*argp_mean)-0.5*(1-5*incl_B)*sin(2*u_mean)+(7/6*incl_B ...
            -1/6*e_mean^2*(1-19/4*incl_B))/e_mean*sin(3*nu_mean+2*argp_mean) ...
            +0.75*incl_B*sin(4*nu_mean+2*argp_mean)+e_mean/8*incl_B*sin(5*nu_mean+2*argp_mean));
        argp_mean = argp_osc - dargp_sp;

        % --- mean anomaly correction ---
        dM_sp = 1.5*j2*(req/p_mean)^2*sqrt(1-e_mean^2)/e_mean*(-(1-3*incl_B)*((1-e_mean^2/4) ...
            *sin(nu_mean)+e_mean/2*sin(2*nu_mean)+e_mean^2/12*sin(3*nu_mean)) ...
            +incl_B*(0.5*(1+1.25*e_mean^2)*sin(nu_mean+2*argp_mean)-e_mean^2/8*sin(nu_mean-2*argp_mean) ...
            -7/6*(1-e_mean^2/28)*sin(3*nu_mean+2*argp_mean)-0.75*e_mean*sin(4*nu_mean+2*argp_mean) ...
            -e_mean^2/8*sin(5*nu_mean+2*argp_mean)));
        M_mean = M_osc - dM_sp;

        % Re-derive true anomaly and argument of latitude from the
        % refined mean anomaly, ready for the next iteration.
        [~, nu_mean] = solve_kepler_equation(M_mean, e_mean);
        u_mean = mod(argp_mean + nu_mean, two_pi);
    end
    % end of the 5-pass fixed-point loop
end

%% ------------------------------------------------------------------
%  Step 3 — assemble the output
%  ------------------------------------------------------------------
mean_elements = zeros(7, 1);
mean_elements(1) = a_mean;
mean_elements(2) = e_mean;
mean_elements(3) = i_mean;
mean_elements(4) = mod(argp_mean, two_pi);
mean_elements(5) = mod(raan_mean, two_pi);
mean_elements(6) = mod(M_mean, two_pi);

% Also report the mean true anomaly for convenience, derived from the
% mean anomaly via Kepler's equation.
[~, nu_out] = solve_kepler_equation(mean_elements(6), e_mean);
mean_elements(7) = nu_out;

end


%% ==========================================================================
function angle = four_quadrant_atan(sin_val, cos_val)
% FOUR_QUADRANT_ATAN  Four-quadrant inverse tangent that always returns
%        an angle in [0, 2*pi), given the sine and cosine of that angle.
%        (Standard atan2 returns values in (-pi, pi]; this remaps to the
%        [0, 2*pi) convention used throughout osc2mean.)
%
% input:  sin_val = sine of angle, cos_val = cosine of angle
% output: angle in radians, 0 <= angle <= 2*pi

epsilon  = 0.0000000001;
half_pi  = 0.5 * pi;

if (abs(sin_val) < epsilon)
    angle = (1 - sign(cos_val)) * half_pi;
    return;
else
    quadrant_base = (2 - sign(sin_val)) * half_pi;
end

if (abs(cos_val) < epsilon)
    angle = quadrant_base;
    return;
else
    angle = quadrant_base + sign(sin_val) * sign(cos_val) * (abs(atan(sin_val / cos_val)) - half_pi);
end

end


%% ==========================================================================
function [ecc_anom, true_anom] = solve_kepler_equation(mean_anom, ecc)
% SOLVE_KEPLER_EQUATION  Solve Kepler's equation (M = E - e*sin(E)) for
%        the eccentric anomaly, then convert to true anomaly. Handles
%        circular, elliptic, and hyperbolic orbits. Uses Danby's method
%        — a higher-order (cubic-convergent) refinement of Newton's
%        method, which converges in very few iterations.
%
% input:   mean_anom = mean anomaly [rad], ecc = eccentricity [-]
% output:  ecc_anom = eccentric anomaly [rad], true_anom = true anomaly [rad]

tol    = 1.0e-10;   % convergence tolerance
two_pi = 2.0 * pi;

M_wrapped = mean_anom - two_pi * fix(mean_anom/two_pi);   % wrap into [0, 2*pi)

% ---- initial guess for the eccentric anomaly ----
if (ecc == 0.0)
    % circular orbit: eccentric anomaly = true anomaly = mean anomaly
    true_anom = M_wrapped;
    ecc_anom  = M_wrapped;
    return;
elseif (ecc < 1.0)
    % elliptic orbit: standard starting guess
    ecc_anom = M_wrapped + 0.85 * sign(sin(M_wrapped)) * ecc;
else
    % hyperbolic orbit: standard starting guess
    ecc_anom = log(2.0 * M_wrapped / ecc + 1.8);
end

% ---- Danby's iteration ----
iter_count = 0;
while (1)

    if (ecc < 1)
        % elliptic orbit: residual(E) = E - e*sin(E) - M = 0
        ecc_sinE = ecc * sin(ecc_anom);
        ecc_cosE = ecc * cos(ecc_anom);
        residual     = ecc_anom - ecc_sinE - M_wrapped;
        d_residual   = 1 - ecc_cosE;
        dd_residual  = ecc_sinE;
        ddd_residual = ecc_cosE;
    else
        % hyperbolic orbit: residual(E) = e*sinh(E) - E - M = 0
        ecc_sinhE = ecc * sinh(ecc_anom);
        ecc_coshE = ecc * cosh(ecc_anom);
        residual     = ecc_sinhE - ecc_anom - M_wrapped;
        d_residual   = ecc_coshE - 1;
        dd_residual  = ecc_sinhE;
        ddd_residual = ecc_coshE;
    end

    iter_count = iter_count + 1;

    if (abs(residual) <= tol || iter_count > 20)
        break;
    end

    % Danby's cubic-convergence correction (successively refined step)
    step1 = -residual / d_residual;
    step2 = -residual / (d_residual + 0.5 * step1 * dd_residual);
    step3 = -residual / (d_residual + 0.5 * step2 * dd_residual ...
             + step2 * step2 * ddd_residual / 6);

    ecc_anom = ecc_anom + step3;
end

if (iter_count > 20)
    fprintf('\n\n   more than 20 iterations in solve_kepler_equation \n\n');
end

% ---- eccentric anomaly -> true anomaly ----
if (ecc < 1)
    sin_true = sqrt(1 - ecc * ecc) * sin(ecc_anom);
    cos_true = cos(ecc_anom) - ecc;
else
    sin_true = sqrt(ecc * ecc - 1) * sinh(ecc_anom);
    cos_true = ecc - cosh(ecc_anom);
end

true_anom = four_quadrant_atan(sin_true, cos_true);

end
