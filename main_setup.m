clear
clc
%constants
G  = 6.6743 * 10^(-11); %N.m^2/kg^2
earth_mass = 5.972 * 10^24; %kg
mu = G*earth_mass; %m^3/s^2
J_2     =  1.08263 * 10^(-3);
%constant for now
mean_rho_500km = 7.30*10^(-13); %kg/m^3
rho_follower = mean_rho_500km;
rho_leader   = mean_rho_500km;

m_f = 2;
m_l = 2;
A_l_min = 0.01; %m^2
A_f_min = 0.01; %m^2
A_l_max = 0.060827625302978; %m^2
A_f_max = 0.060827625302978; %m^2

B_f_min = 0.013012490810596; %(C_d*A_f_min)/m_f;
B_f_max = 0.078628595978917; %(C_d*A_f_max)/m_f;

B_l_min = 0.013012490810596; %(C_d*A_l_min)/m_l;
B_l_max = 0.078628595978917; %(C_d*A_l_max)/m_l;


%initial conditions Orbit
earth_radi      = 6371*10^3; %m

a_leader        = 6786 *10^3; %m 
a_follower      = 6786 *10^3; %m

e_leader       = 0.01; %eccenticity
e_follower     = 0.01;  %eccentricity

i_leader        = 20; % inclination
i_follower      = 20; % inclination

Omega_leader    =  95; % right ascension of accending node
Omega_follower  =  95; % right ascension of accendng node

omega_leader    = 25; % argument of periapsis;
omega_follower  = 25; % argument of periapsis;

TA_leader     = 100.000;   %True anomaly
TA_follower   = 100.001;  %True anomaly
%% -----------


n_leader    = sqrt(mu/a_leader^3);   %mean angular motion
n_follower  = sqrt(mu/a_follower^3); %mean angular motion

Period_leader = 2*pi*sqrt((a_leader)^3/mu);
Period_follower = 2*pi*sqrt((a_follower)^3/mu);

%% Control Refrence
delta_a_ref = 0;
delta_lambda_ref = 20/a_leader;