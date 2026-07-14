function yhat = skewHistModel_A(p, x)
% skewHistModel_A
% Histogram-normalized skew model (Gaussian + one-sided exponential tail) + baseline.
%
% This model is designed for fitting *counts/bin* data directly, without
% multiplying by bin widths. The "shape" is computed at bin centers and then
% normalized so that sum(shape) = 1 over the provided bins.
%
% Then:
%   yhat(i) = b0 + A * w(i)
% where w is dimensionless and sums to 1, so A is the total counts above
% baseline assigned to the peak across the bins included in the fit.
%
% p = [A, mu, sigma, b0, f, tau]
% x = bin start positions (monotonically increasing; uniform or non-uniform ok)
%
% Output:
%   yhat = predicted counts/bin at each bin (same size as x)

    % ---- unpack ----
    A     = p(1);
    mu    = p(2);
    sigma = abs(p(3)) + eps;
    b0    = p(4);
    f     = p(5);
    tau   = abs(p(6)) + eps;

    % ---- sanitize inputs ----
    x = x(:);
    if numel(x) < 3
        yhat = nan(size(x));
        return;
    end

    % Ensure monotonic increasing
    dx = diff(x);
    if any(~isfinite(dx)) || any(dx <= 0)
        error('skewHistModel_A:AxisNotMonotonic', ...
              'x must be monotonically increasing bin start positions.');
    end

    % Bin centers (use last dx for last bin like your original)
    dx = [dx; dx(end)];
    xc = x + 0.5*dx;

    % ---- unnormalized "shape" at bin centers ----
    % Gaussian (pdf units 1/x)
    g = (1 ./ (sigma * sqrt(2*pi))) .* exp(-0.5 .* ((xc - mu) ./ sigma).^2);

    % One-sided exponential tail starting at mu (also 1/x)
    u = xc - mu;
    e = (1 ./ tau) .* exp(-max(u,0) ./ tau) .* (u >= 0);

    % Mixture (still 1/x)
    s = (1 - f) .* g + f .* e;

    % ---- histogram normalization across bins ----
    % Convert to dimensionless weights that sum to 1.
    % IMPORTANT: we do NOT multiply by dx here (by design).
    ssum = sum(s);
    if ~isfinite(ssum) || ssum <= 0
        % Degenerate parameter region -> return baseline only
        yhat = b0 * ones(size(x));
        return;
    end
    w = s ./ ssum;   % dimensionless, sum(w)=1

    % ---- model in counts/bin ----
    yhat = b0 + A .* w;
end