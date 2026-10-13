
%% Controller Performance Evaluation

% Simulation time
t = out.tout;

% Relative orbital elements
delta_a = out.states.Data(:,1);
delta_lambda = out.states.Data(:,2);

% Ballistic coefficients
beta_l = out.ballistic.Data(:,4);
beta_f = out.ballistic.Data(:,3);

% Desired formation
delta_a_ref = 0;
delta_lambda_ref = 20/a_leader;

% Ballistic coefficient constraints
beta_min = B_f_min;
beta_max = B_f_max;

% Settling tolerances
eps_a = 0.5/a_leader; %meter
eps_lambda = 2/a_leader;  %meter

% Calculate performance metrics
metrics = evaluate_orbit_control( ...
    t, a_leader,delta_a, delta_lambda, ...
    beta_l, beta_f, ...
    delta_a_ref, delta_lambda_ref, ...
    beta_min, beta_max, ...
    eps_a, eps_lambda);

%% Display results

fprintf('\n========== CONTROL PERFORMANCE ==========\n');

fprintf('Normalized Effort (Leader): %.6f\n', ...
    metrics.J_beta_N_leader);

fprintf('Normalized Effort (Follower): %.6f\n', ...
    metrics.J_beta_N_follower);

fprintf('Normalized Effort (Mean): %.6f\n', ...
    metrics.J_beta_N_mean);

fprintf('Settling Time delta_a: %.2f Orbit\n', ...
    metrics.Ts_delta_a/Period_leader);

fprintf('Settling Time delta_lambda: %.2f Orbit\n', ...
    metrics.Ts_delta_lambda/Period_leader);

fprintf('Maximum Overshoot delta_a: %.2f %%\n', ...
    metrics.Mp_delta_a_percent);

fprintf('Maximum Overshoot delta_lambda: %.2f %%\n', ...
    metrics.Mp_delta_lambda_percent);

fprintf('RMSE a*delta_a: %.3f m\n', ...
    metrics.rmse_delta_a*a_leader);

fprintf('RMSE a*delta_lambda: %.3f m\n', ...
    metrics.rmse_delta_lambda*a_leader);

fprintf('Final Error a*delta_a: %.3f m\n', ...
    metrics.final_e_delta_a*a_leader);

fprintf('Final Error a*delta_lambda: %.3f m\n', ...
    metrics.final_e_delta_lambda*a_leader);
