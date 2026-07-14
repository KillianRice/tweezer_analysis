function chat = skewRateModel_A(p, x)
% Normalized Gaussian + normalized one-sided exponential tail + baseline
% Works with non-uniform bin widths (e.g., voltage from exponential ramp)
%
% p = [A, mu, sigma, b0, f, tau]
% x = bin start positions (must be monotonically increasing, not necessarily uniform)

    A     = p(1);
    mu    = p(2);
    sigma = abs(p(3)) + eps;
    b0    = p(4);
    f     = p(5);
    tau   = abs(p(6)) + eps;

    x = x(:);

    % Per-bin widths (non-uniform allowed)
    dx = diff(x);
    if any(~isfinite(dx)) || any(dx <= 0)
        error('Axis must be monotonically increasing.');
    end
    dx = [dx; dx(end)];          % last bin: reuse previous width

    % Bin centers
    xc = x + 0.5*dx;

    % ---- Normalized Gaussian PDF (area = 1) ----
    g = (1 ./ (sigma * sqrt(2*pi))) .* exp(-0.5 .* ((xc - mu) ./ sigma).^2);

    % ---- Normalized one-sided exponential PDF (area = 1) ----
    u = xc - mu;
    e = (1 ./ tau) .* exp(-max(u,0) ./ tau) .* (u >= 0);

    % ---- Mixture PDF (still normalized to 1) ----
    pdf = (1 - f) .* g + f .* e;

    % ---- Convert PDF (1/x-unit) to counts per bin ----
    chat = b0 + A .* (pdf .* dx);
end