function [V, dV] = timeToVoltageAxis(t)
% timeToVoltageAxis
% Converts monotonically increasing time bins to voltage bins
% using the calibrated exponential ramp:
%
%   V(t) = Vinf * (1 - exp(-t/tauRC))
%
% Voltage is returned as positive magnitude starting at 0,
% scaled to mV.
%
% Outputs:
%   V  : voltage bin start positions (mV)
%   dV : per-bin voltage widths (mV)

    % -------------------------------
    % HARD-CODED RAMP CALIBRATION
    % -------------------------------
    Vinf  = -0.140912;      % asymptotic voltage (Volts)
    tauRC = 4.81377e-6;     % RC time constant (seconds)
    scale = 1000;           % convert V -> mV

    t = t(:);

    if any(~isfinite(t))
        error('timeToVoltageAxis:NonFiniteTime','t contains non-finite values.');
    end

    % Raw mapping
    V = Vinf .* (1 - exp(-t ./ tauRC));

    % Force positive ramp magnitude starting at 0
    V = abs(V - V(1));

    % Convert to mV
    V = V .* scale;

    % Compute per-bin widths (non-uniform)
    dV = diff(V);

    if any(~isfinite(dV)) || any(dV <= 0)
        error('timeToVoltageAxis:NonMonotonicV', ...
              'Voltage axis not strictly increasing after mapping.');
    end

    % Repeat last bin width
    dV = [dV; dV(end)];
end