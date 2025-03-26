function [modDropFit] = DropTimeCalibration(time,pos)
% Function to fit drop time absorption or fluorescence data.
%
% INPUTS:
%   time - Drop time in seconds
%   pos  - spatial position (arb. unit)

%% Fitting function
s2 = @(coeffs,t) (coeffs(2) - 0.49.*coeffs(1).*(t.^2)); % Y0 - C*(1/2gt^2), where C is in pixels/meter.

%% Guessing
%[peakTime,peakLoc] = findpeaks(smooth(pos),'minpeakdistance',2,'minpeakheight',mean(pos));

offsetGuess = max(pos);
CGuess = 10^5; %pixels/meter

%% Fitting
%statOpt = statset('RobustWgtFun','bisquare');
[modDropFit] = nlinfit(time,pos,s2,[CGuess offsetGuess]);

%% Print and Plot results
fprintf('\n Top pixel %g px\n\n',modDropFit(2))
fprintf('\n Camera pixel calibration is %g m/px\n\n',1/modDropFit(1))
timeInterp = linspace(time(1),time(end),500);
dataH = plot(time,pos,'ro',timeInterp,s2(modDropFit,timeInterp),'b--');
grid on; axis tight
xlabel('Time [s]','FontSize',24,'FontWeight','Bold')
ylabel('Center Position','FontSize',24,'FontWeight','Bold')
set(dataH(1),'MarkerSize',6,'MarkerFaceColor','r')
set(dataH(2),'LineWidth',2)
end