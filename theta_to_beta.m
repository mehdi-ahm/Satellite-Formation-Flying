function [beta,C_d,A] = theta_to_beta(theta)
% THETA_TO_BETA  Convert a satellite pitch angle theta [rad]
% to the corresponding ballistic coefficient beta [m^2/kg].
%
% 
%
%   This is a direct forward evaluation (no lookup/interpolation needed,
%   since beta(theta) is single-valued) using the same flat-plate /
%   area model as beta_to_theta's table build.

    m = 2; % spacecraft mass [kg]

    % Clamp theta
    theta_min = 0;
    theta_max = 83.6557*pi/180;   % [rad], peak of beta(theta)

    theta = min(max(theta, theta_min), theta_max);
    % C_d model constant (NRLMSISE-00 based)
    alpha = (7.5*10^(-17)*5.14*10^7*10^6*1011) / (1.0 + (7.5*10^(-17)*5.14*10^7*10^6*1011));
    r = sqrt(1-alpha);

    %theta pitch of satellite [rad]

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

    A = A_ref_solar + A_ref_box_top + A_ref_box_front;

    % Total drag coefficient (area-weighted blend)
    C_d = (A_ref_box_front/A)*C_d_flat_theta_inv ...
              + ((A-A_ref_box_front)/A)*C_d_flat_theta;

    beta = (A * C_d) / m;

end