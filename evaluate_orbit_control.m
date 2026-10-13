
function metrics = evaluate_orbit_control(t,a_leader, delta_a, delta_lambda, beta_l, beta_f, ...
    delta_a_ref, delta_lambda_ref, beta_min, beta_max, eps_a, eps_lambda)
%EVALUATE_ORBIT_CONTROL Quantitative performance metrics for formation control.
%
% INPUTS
%   t                 N-by-1 simulation time [s], strictly increasing
%   delta_a           N-by-1 relative semi-major axis [-]
%   delta_lambda      N-by-1 relative mean longitude [rad]
%   beta_l, beta_f    N-by-1 ballistic coefficients [m^2/kg]
%   delta_a_ref       scalar or N-by-1 reference [-]
%   delta_lambda_ref  scalar or N-by-1 reference [rad]
%   beta_min/max      ballistic coefficient limits [m^2/kg]
%   eps_a             settling tolerance for delta_a [-]
%   eps_lambda        settling tolerance for delta_lambda [rad]
%
% OUTPUT
%   metrics: structure containing performance metrics.
%
% NORMALIZED CONTROL EFFORT
%   J_beta,N = (1/T) * integral((beta-beta_min)/(beta_max-beta_min) dt)
%
% SETTLING TIME
%   Time after which the tracking error remains within tolerance.
%
% MAXIMUM OVERSHOOT
%   Percentage overshoot relative to a nonzero constant reference.
%   Returns NaN for zero or time-varying references.
%
% NOTE:
%   delta_lambda uses direct subtraction, assuming a continuous
%   (unwrapped) angular branch.

    %% 1. Prepare inputs

    t = t(:);
    delta_a = delta_a(:);
    delta_lambda = delta_lambda(:);
    beta_l = beta_l(:);
    beta_f = beta_f(:);

    n = numel(t);

    if n < 2 || any(~isfinite(t)) || any(diff(t) <= 0)
        error('t must contain at least two finite, strictly increasing samples.');
    end

    series = {delta_a, delta_lambda, beta_l, beta_f};

    if any(cellfun(@(x) numel(x) ~= n || any(~isfinite(x)), series))
        error('All simulated signals must be finite and the same length as t.');
    end

    if ~isscalar(beta_min) || ~isscalar(beta_max) || ...
            ~isfinite(beta_min) || ~isfinite(beta_max) || beta_max <= beta_min
        error('Require finite beta_max > beta_min.');
    end

    if ~isscalar(eps_a) || ~isscalar(eps_lambda) || ...
            ~isfinite(eps_a) || ~isfinite(eps_lambda) || eps_a < 0 || eps_lambda < 0
        error('Settling tolerances must be finite nonnegative scalars.');
    end

    delta_a_ref = expandReference(delta_a_ref, n, 'delta_a_ref');
    delta_lambda_ref = expandReference(delta_lambda_ref, n, 'delta_lambda_ref');


    %% 2. Relative orbital element errors

    metrics.e_delta_a = delta_a - delta_a_ref;

    metrics.e_aXdelta_a = a_leader*(delta_a - delta_a_ref);
    
    metrics.e_delta_lambda = delta_lambda - delta_lambda_ref;

    metrics.e_aXdelta_lambda = a_leader*(delta_lambda - delta_lambda_ref);

    metrics.abs_e_delta_a = abs(metrics.e_delta_a);

    metrics.abs_e_delta_lambda = abs(metrics.e_delta_lambda);


    %% 3. Root Mean Square Error (RMSE)

    T = t(end) - t(1);

    metrics.rmse_delta_a = sqrt( ...
        trapz(t, metrics.e_delta_a.^2) / T);

    metrics.rmse_delta_lambda = sqrt( ...
        trapz(t, metrics.e_delta_lambda.^2) / T);


    %% 4. Normalized Control Effort

    beta_range = beta_max - beta_min;

    metrics.u_normalized_leader = ...
        (beta_l - beta_min) / beta_range;

    metrics.u_normalized_follower = ...
        (beta_f - beta_min) / beta_range;

    metrics.J_beta_N_leader = ...
        trapz(t, metrics.u_normalized_leader) / T;

    metrics.J_beta_N_follower = ...
        trapz(t, metrics.u_normalized_follower) / T;

    metrics.J_beta_N_mean = ...
        0.5 * (metrics.J_beta_N_leader + metrics.J_beta_N_follower);

    % Check admissible ballistic coefficient limits

    metrics.beta_bounds_respected = ...
        all(beta_l >= beta_min & beta_l <= beta_max) && ...
        all(beta_f >= beta_min & beta_f <= beta_max);


    %% 5. Settling Time

    metrics.Ts_delta_a = settlingTime( ...
        t, metrics.abs_e_delta_a, eps_a);

    metrics.Ts_delta_lambda = settlingTime( ...
        t, metrics.abs_e_delta_lambda, eps_lambda);


    %% 6. Maximum Overshoot (%)

    metrics.Mp_delta_a_percent = ...
        overshoot(delta_a, delta_a_ref);

    metrics.Mp_delta_lambda_percent = ...
        overshoot(delta_lambda, delta_lambda_ref);


    %% 7. Final Tracking Errors

    metrics.final_e_delta_a = metrics.e_delta_a(end);

    metrics.final_e_delta_lambda = metrics.e_delta_lambda(end);

    metrics.settled_delta_a = isfinite(metrics.Ts_delta_a);

    metrics.settled_delta_lambda = isfinite(metrics.Ts_delta_lambda);

end


%% ========================================================================
%  LOCAL FUNCTION: REFERENCE HANDLING
%  ========================================================================

function ref = expandReference(ref, n, name)

    if isscalar(ref)
        ref = repmat(ref, n, 1);
    else
        ref = ref(:);
    end

    if numel(ref) ~= n || any(~isfinite(ref))
        error('%s must be a finite scalar or have N samples.', name);
    end

end


%% ========================================================================
%  LOCAL FUNCTION: SETTLING TIME
%  ========================================================================

function Ts = settlingTime(t, absError, tolerance)

    idx = find(absError > tolerance, 1, 'last');

    if isempty(idx)

        % Already inside tolerance throughout simulation
        Ts = 0;

    elseif idx == numel(t)

        % Never settled during simulation
        Ts = NaN;

    else

        % First sample after the last tolerance violation
        Ts = t(idx+1) - t(1);

    end

end


%% ========================================================================
%  LOCAL FUNCTION: MAXIMUM OVERSHOOT
%  ========================================================================

function mp = overshoot(y, ref)

    r = ref(end);

    % Overshoot percentage undefined for zero or varying reference
    if abs(r) <= 1e-12 || ...
            any(abs(ref-r) > 1e-10*max(1,abs(r)))

        mp = NaN;
        return;

    end

    initialError = y(1) - r;

    if abs(initialError) <= 1e-12

        mp = NaN;

    elseif initialError < 0

        % Initial output is below reference
        mp = 100 * max(0, max(y-r)) / abs(r);

    else

        % Initial output is above reference
        mp = 100 * max(0, max(r-y)) / abs(r);

    end

end
