function sigma_est = estimateSigmaLeftHalfMax(t, counts)
% Estimate Gaussian sigma from the LEFT half-maximum width.
%
% Uses only the left side of the peak to avoid distortion
% from a long right-side exponential tail.
%
% Inputs:
%   t       : time vector (bin start times, column vector)
%   counts  : counts per bin (same length as t)
%
% Output:
%   sigma_est : estimated Gaussian sigma (seconds)

    % Ensure column vectors
    t = t(:);
    counts = counts(:);

    N = numel(t);
    if N < 3
        error('Not enough data points to estimate sigma.');
    end

    % --- Estimate baseline from early bins ---
    nb = min(20, N);  % use first 20 bins (or fewer if short trace)
    b0_est = median(counts(1:nb));

    % --- Find peak ---
    [peakVal, ip] = max(counts);
    mu_est = t(ip);

    % --- Compute half-maximum level above baseline ---
    yhalf = b0_est + 0.5 * (peakVal - b0_est);

    % --- Find last point on left below half-max ---
    leftIdx = find(counts(1:ip) <= yhalf, 1, 'last');

    if isempty(leftIdx) || leftIdx >= ip
        % Fallback if half-max crossing not found
        dt = median(diff(t));
        sigma_est = 3 * dt;
        return
    end

    % --- Linear interpolation for better half-max crossing estimate ---
    t1 = t(leftIdx);
    t2 = t(leftIdx + 1);
    y1 = counts(leftIdx);
    y2 = counts(leftIdx + 1);

    % Protect against division issues
    if y2 == y1
        dt = median(diff(t));
        sigma_est = 3 * dt;
        return
    end

    % Interpolated crossing time
    t_half = t1 + (yhalf - y1) * (t2 - t1) / (y2 - y1);

    % Left half-width
    dt_left = mu_est - t_half;

    % Convert half-width to sigma for Gaussian:
    % FWHM = 2*sqrt(2 ln 2)*sigma
    % Here we only have half-width, so:
    sigma_est = dt_left / sqrt(2*log(2));

    % Enforce minimum sigma ~ one bin width
    dt = median(diff(t));
    sigma_est = max(sigma_est, dt);

end