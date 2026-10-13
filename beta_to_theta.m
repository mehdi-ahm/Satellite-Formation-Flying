function theta_rad = beta_to_theta(beta_c)
% BETA_TO_THETA  Convert a ballistic coefficient to a pitch angle (rad)
%                using the rising branch of the beta(theta) curve.
%
%   theta_deg = beta_to_theta(beta_c)
%
%   The beta(theta) lookup table is built ONCE (on the first call) and
%   cached in persistent variables. Every subsequent call only does a
%   1-D interpolation, so this is cheap to call repeatedly (e.g. inside
%   a Simulink MATLAB Function block or NMPC cost/constraint function).

    persistent beta_rising theta_rising B_min B_max initialized

    if isempty(initialized)
        [beta_rising, theta_rising, B_min, B_max] = build_beta_theta_table();
        initialized = true;
    end

    if beta_c < B_min || beta_c > B_max
        % Clamp instead of erroring - keeps this safe for Simulink
        % compile-time probing (zero/default inputs) and NMPC solver
        % excursions outside the achievable range.
        beta_c = min(max(beta_c, B_min), B_max);
    end

    theta_rad = interp1(beta_rising, theta_rising, beta_c, 'linear');
    theta_deg = theta_rad * 180/pi;

end


function [beta_rising, theta_rising, B_min, B_max] = build_beta_theta_table()

    n = 4; % power number of data
    num_points = 90*10^n + 1;

    m = 2; % spacecraft mass [kg]

    % C_d model constant (NRLMSISE-00 based)
    alpha = (7.5*10^(-17)*5.14*10^7*10^6*1011) / (1.0 + (7.5*10^(-17)*5.14*10^7*10^6*1011));
    r = sqrt(1-alpha);

    C_d_vec       = zeros(num_points,1);
    theta_vec     = zeros(num_points,1);
    A_ref_tot_vec = zeros(num_points,1);

    for k = 1:num_points

        theta = (pi/(180*10^n))*(k-1); % pitch of satellite [rad]

        % flat plate at incidence theta
        C_d_flat_theta     = 2*(1+(2/3)*r*sin(theta));
        C_d_flat_theta_inv = 2*(1+(2/3)*r*sin(pi/2-theta));

        % Areas
        % 2 solar panels
        A_tot_solar = 2*15*15/10^4;
        A_ref_solar = A_tot_solar*sin(theta);

        % box (front and upper side)
        A_tot_box_front = 10*10/10^4;
        A_tot_box_top   = 15*10/10^4;

        A_ref_box_front = A_tot_box_front*cos(theta);
        A_ref_box_top   = A_tot_box_top*sin(theta);

        A_ref_tot = A_ref_solar + A_ref_box_top + A_ref_box_front;

        % Total drag coefficient (area-weighted blend)
        C_d_total = (A_ref_box_front/A_ref_tot)*C_d_flat_theta_inv ...
                  + ((A_ref_tot-A_ref_box_front)/A_ref_tot)*C_d_flat_theta;

        C_d_vec(k)       = C_d_total;
        theta_vec(k)     = theta;
        A_ref_tot_vec(k) = A_ref_tot;
    end

    beta_vec = (A_ref_tot_vec.*C_d_vec)/m;

    % Restrict to the rising branch only: theta = 0 to peak
    [beta_pk, idx_pk] = max(beta_vec);

    theta_rising = theta_vec(1:idx_pk);
    beta_rising  = beta_vec(1:idx_pk);

    B_min = min(beta_rising);
    B_max = beta_pk;

end