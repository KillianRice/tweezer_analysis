function [modDropFit] = DropTimeCalibrationErr(time,pos,avgPlot,err)
% Function to fit drop time absorption or fluorescence data.
%
% INPUTS:
%   time - Drop time in seconds
%   pos  - spatial position (arb. unit)

%% Fitting function
s2 = @(coeffs,t) (coeffs(2) + coeffs(3).*t - 0.5*9.81*(0.926)*coeffs(1).*(t.^2)); % Y0 - C*(1/2gt^2) + Vt, where C is in pixels/meter.


t = time;
y = pos;
sigma = err;

%% Initial guesses
beta0 = [1e5, max(y), 0];

%% Weighted fitting
residual = @(beta) (y - s2(beta,t)) ./ sigma;

opts = optimoptions('lsqnonlin','Display','off');
[modDropFit, resnorm, residuals, exitflag, output, lambda, J] = ...
    lsqnonlin(residual, beta0, [], [], opts);

%% Print and Plot results
if avgPlot
    fprintf('Avg Plot Stats')
end
fprintf('\n Top pixel %g px\n\n',modDropFit(2))
fprintf('\n Camera pixel calibration is %g m/px\n\n',1/modDropFit(1))
fprintf('\n Sample Velocity %g px/s\n\n',modDropFit(3))
timeInterp = linspace(time(1),time(end),500);
dataH = plot(time,pos,'ro',timeInterp,s2(modDropFit,timeInterp),'b--');
grid on; axis tight
if avgPlot
    title('Average Center Position vs Time')
end
xlabel('Time [s]','FontSize',24,'FontWeight','Bold')
ylabel('Center Position','FontSize',24,'FontWeight','Bold')
set(dataH(1),'MarkerSize',6,'MarkerFaceColor','r')
set(dataH(2),'LineWidth',2)

%% Plot with err
% Fitted coefficients
Cfit      = 1/modDropFit(1);
offsetFit = modDropFit(2);
velFit    = modDropFit(3);

errorbar(time, pos, err, 'o');
hold on

tFit = linspace(min(time), max(time), 200);
plot(tFit, s2(modDropFit,tFit), 'r-', 'LineWidth', 2);
title('Average Center Position vs Time')
xlabel('Time (s)');
ylabel('Position (Px)');
legend('Data \pm SEM', 'Weighted fit');
eqnText = sprintf(['y(t) = A + V*t + 0.5g*C*t^2\n' ...
                   'A (px) = %.4g\n' 'V (px/s) = %.4g\n' '1/C (m/px) = %.4g'], ...
                   offsetFit, velFit, Cfit);

text(0.05, 0.95, eqnText, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 11, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');

end