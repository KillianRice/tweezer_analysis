function [modDropFit] = DropTimeCalibration(time,pos,avgPlot)
% Function to fit drop time absorption or fluorescence data.
%
% INPUTS:
%   time - Drop time in seconds
%   pos  - spatial position (arb. unit)

%% Fitting function
s2 = @(coeffs,t) (coeffs(2) + coeffs(3).*t - 4.9.*coeffs(1).*(t.^2)); % Y0 - C*(1/2gt^2) + Vt, where C is in pixels/meter.

%% Guessing
%[peakTime,peakLoc] = findpeaks(smooth(pos),'minpeakdistance',2,'minpeakheight',mean(pos));

offsetGuess = max(pos);
CGuess = 10^5; %pixels/meter
vel = 0;

%% Fitting
%statOpt = statset('RobustWgtFun','bisquare');
[modDropFit] = nlinfit(time,pos,s2,[CGuess offsetGuess vel]);

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
end