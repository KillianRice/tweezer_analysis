function funcOut = fit_release_recapture_AveragedCounts( ...
analyVar, indivDataset, avgDataset, useroptions)
%FIT_RELEASE_RECAPTURE_AVERAGEDCOUNTS
% Fit release-and-recapture data using averaged raw counts.
%
% Workflow:
% 1. Average the data for each tweezer ROI.
% 2. Fit each scan using a Monte Carlo recapture model.
% 3. Combine corresponding points from all tweezers.
% 4. Fit the combined all-tweezer data.
% 5. Plot the individual and/or combined fits.
%
% Model:
%
% Counts(t) = Background + Amplitude * P_recapture(t,T)
%
% The Monte Carlo samples are generated once and reused for every fit.
% This keeps the Monte Carlo noise fixed while the fit changes T.
%
% See create_default_options() below for available options.

%% Options

if nargin < 4
    useroptions = struct;
end

options = create_default_options();

% User-supplied options overwrite the defaults.
names = fieldnames(useroptions);
for i = 1:numel(names)
    options.(names{i}) = useroptions.(names{i});
end

if isempty(options.XAxisLabel)
    options.XAxisLabel = options.IndVarField;
end

if isempty(options.YAxisLabel)
    options.YAxisLabel = options.DepVarField;
end

%% Basic setup

scanIDs = analyVar.uniqScanList(:);
numScanIDs = numel(scanIDs);

if analyVar.UseTweezer
    roiFile = fullfile(analyVar.analyOutDir,'tweezerROI.mat');
    roiData = load(roiFile,'tweezerROI');
    numTweezers = size(roiData.tweezerROI.centersXY,1);
else
    numTweezers = 1;
end

fprintf('Release-and-recapture fit\n');
fprintf('  Tweezers: %d\n',numTweezers);
fprintf('  Scan IDs: %d\n\n',numScanIDs);

%% Monte Carlo samples

% Use the same random samples for every fit. This is important:
% otherwise the objective function itself would change randomly
% during optimization.
oldRNG = rng;
cleanupRNG = onCleanup(@() rng(oldRNG)); %#ok<NASGU>

rng(options.RandomSeed,'twister');

positionSamples = randn(options.NumParticles,3);
velocitySamples = randn(options.NumParticles,3);

physicalModel = create_physical_model(options);

%% Storage

xByTweezer = cell(numTweezers,numScanIDs);
yByTweezer = cell(numTweezers,numScanIDs);
yerrByTweezer = cell(numTweezers,numScanIDs);

individualFits = cell(numTweezers,numScanIDs);


%% Fit each tweezer
for tweezer = 1:numTweezers

    if analyVar.UseTweezer
        fprintf('Tweezer %d / %d\n',tweezer,numTweezers);
    end

    [xavg,yavg,yerr] = get_averages( ...
        analyVar, ...
        indivDataset, ...
        avgDataset, ...
        options.IndVarField, ...
        options.DepVarField, ...
        options.Statistics, ...
        tweezer);

    for scan = 1:numScanIDs

        % Clean up this scan and make sure all errors are usable.
        [x,y,error] = clean_scan_data( ...
            xavg{scan}, ...
            yavg{scan}, ...
            yerr{scan}, ...
            options.MinimumError);

        xByTweezer{tweezer,scan} = x;
        yByTweezer{tweezer,scan} = y;
        yerrByTweezer{tweezer,scan} = error;

        % Not enough points for a meaningful 3-parameter fit.
        if numel(x) < 3
            continue;
        end
        
        if ~analyVar.SkipIndivTwzrFitting
            disp("Doing Indiv Tweezer Fitting")
            fitResult = fit_single_raw_count_scan( ...
                x*options.TimeScaleToSeconds, ...
                y, ...
                error, ...
                positionSamples, ...
                velocitySamples, ...
                physicalModel, ...
                options);
    
            fitResult.xDataOriginalUnits = x;
            fitResult.scanID = scanIDs(scan);
            fitResult.tweezerNum = tweezer;
    
            individualFits{tweezer,scan} = fitResult;
    
            if options.PlotIndividualTweezers
                plot_single_fit( ...
                    x,y,error,fitResult, ...
                    scanIDs(scan),tweezer,options,false);
            end
        end
    end
end

%% Combine corresponding points from all tweezers

combinedX = cell(numScanIDs,1);
combinedY = cell(numScanIDs,1);
combinedYerr = cell(numScanIDs,1);
combinedFits = cell(numScanIDs,1);

for scan = 1:numScanIDs

    [combinedX{scan}, ...
     combinedY{scan}, ...
     combinedYerr{scan}] = combine_tweezer_scans( ...
        xByTweezer(:,scan), ...
        yByTweezer(:,scan), ...
        yerrByTweezer(:,scan));

    if numel(combinedX{scan}) < 3
        continue;
    end

    fitResult = fit_single_raw_count_scan( ...
        combinedX{scan}*options.TimeScaleToSeconds, ...
        combinedY{scan}, ...
        combinedYerr{scan}, ...
        positionSamples, ...
        velocitySamples, ...
        physicalModel, ...
        options);

    fitResult.xDataOriginalUnits = combinedX{scan};
    fitResult.scanID = scanIDs(scan);
    fitResult.tweezerNum = NaN;

    combinedFits{scan} = fitResult;

    if options.PlotAllTweezerAverage
        plot_single_fit( ...
            combinedX{scan}, ...
            combinedY{scan}, ...
            combinedYerr{scan}, ...
            fitResult, ...
            scanIDs(scan), ...
            NaN, ...
            options, ...
            true);
    end
end

%% Optional overview plot

if options.PlotAllFitsTogether
    plot_all_combined_fits( ...
        combinedX, ...
        combinedY, ...
        combinedYerr, ...
        combinedFits, ...
        scanIDs, ...
        analyVar, ...
        options);
end

%% Return results

funcOut.analyVar = analyVar;
funcOut.indivDataset = indivDataset;
funcOut.avgDataset = avgDataset;

end

%% ========================================================================
% OPTIONS
% ========================================================================

function options = create_default_options()

options.IndVarField = 'imagevcoAtom';
options.DepVarField = 'OD_TotalCounts';
options.Statistics = 'gaussian';

options.XAxisLabel = '';
options.YAxisLabel = '';

options.XAxisScale = 'linear';
options.YAxisScale = 'linear';

% LabVIEW stores the release time in ms.
options.TimeScaleToSeconds = 0.001;

% Plotting
options.PlotIndividualTweezers = false;
options.PlotAllTweezerAverage = true;
options.PlotAllFitsTogether = true;
options.ShowFitAnnotation = true;

options.FitTitle = 'Release-and-Recapture Raw-Count Fit';

% Trap / atom parameters
options.MassAMU = 87.9056;       % Sr-88
options.OmegaXHz = 82e3;
options.OmegaYHz = 82e3;
options.OmegaZHz = 10.5e3;

options.TrapWavelength = 532e-9;

options.WaistX = 1e-6;
options.WaistY = 1e-6;
options.WaistZ = pi*(1e-6)^2/(532e-9);

options.Gravity = 9.8;
options.GravityAxis = 'z';

% Monte Carlo
options.NumParticles = 20000;
options.RandomSeed = 1;

% Fit parameters
options.InitialTemperature = 10e-6;
options.MinimumTemperature = 1e-6;
options.MaximumTemperature = 40e-6;

options.InitialAmplitude = [];
options.InitialBackground = [];

options.MinimumError = 1;
options.NumFitCurvePoints = 50000;

options.FitOptions = optimset( ...
    'Display','off', ...
    'MaxFunEvals',7000, ...
    'MaxIter',10000, ...
    'TolX',1e-7, ...
    'TolFun',1e-7);

end

%% ========================================================================
% PHYSICAL MODEL
% ========================================================================

function physicalModel = create_physical_model(options)

kB = 1.380649e-23;
amu = 1.66053906660e-27;

mass = options.MassAMU*amu;

omegaX = 2*pi*options.OmegaXHz;
omegaY = 2*pi*options.OmegaYHz;
omegaZ = 2*pi*options.OmegaZHz;

k = 2*pi/options.TrapWavelength;

physicalModel.boltzmannConstant = kB;
physicalModel.mass = mass;

physicalModel.omegaX = omegaX;
physicalModel.omegaY = omegaY;
physicalModel.omegaZ = omegaZ;

% Harmonic approximation to the trap depth.
physicalModel.trapDepthX = mass*omegaX^2/(2*k^2);
physicalModel.trapDepthY = mass*omegaY^2/(2*k^2);
physicalModel.trapDepthZ = mass*omegaZ^2/(2*k^2);

end

%% ========================================================================
% DATA CLEANING
% ========================================================================

function [xData,yData,errorData] = ...
clean_scan_data(xData,yData,errorData,minimumError)

xData = xData(:);
yData = yData(:);

if isempty(errorData)
    errorData = nan(size(yData));
else
    errorData = errorData(:);
end

% Remove points with missing x/y values.
valid = isfinite(xData) & isfinite(yData);

xData = xData(valid);
yData = yData(valid);
errorData = errorData(valid);

% Put the release times in increasing order.
[xData,order] = sort(xData);

yData = yData(order);
errorData = errorData(order);

% Replace missing/zero uncertainties with a representative error.
validErrors = errorData(isfinite(errorData) & errorData > 0);

if isempty(validErrors)
    replacementError = minimumError;
else
    replacementError = max(median(validErrors),minimumError);
end

errorData(~isfinite(errorData) | errorData <= 0) = replacementError;
errorData = max(errorData,minimumError);

end

%% ========================================================================
% FIT
% ========================================================================

function fitResult = fit_single_raw_count_scan( ...
releaseTimes, ...
rawCounts, ...
rawCountErrors, ...
positionSamples, ...
velocitySamples, ...
physicalModel, ...
options)

%% Initial guess

if isempty(options.InitialAmplitude)
    initialAmplitude = max(rawCounts)-min(rawCounts);
    initialAmplitude = max(initialAmplitude,1);
else
    initialAmplitude = options.InitialAmplitude;
end

if isempty(options.InitialBackground)
    initialBackground = min(rawCounts);
else
    initialBackground = options.InitialBackground;
end

% Temperature is fitted through a logistic parameter so that it
% automatically stays between the chosen minimum and maximum.
fraction = ...
    (options.InitialTemperature-options.MinimumTemperature) / ...
    (options.MaximumTemperature-options.MinimumTemperature);

fraction = min(max(fraction,1e-9),1-1e-9);

initialTemperatureParameter = log(fraction/(1-fraction));

initialParameters = [ ...
    initialTemperatureParameter
    log(initialAmplitude)
    initialBackground ];

%% Fit form
objective = @(p) raw_count_objective( ...
    p, ...
    releaseTimes, ...
    rawCounts, ...
    rawCountErrors, ...
    positionSamples, ...
    velocitySamples, ...
    physicalModel, ...
    options);

%% Find best fit parameters that minimize the residuals
[bestParameters,objectiveValue,exitFlag,optimizerOutput] = ...
    fminsearch(objective,initialParameters,options.FitOptions);

[temperature,amplitude,background] = ...
    decode_fit_parameters(bestParameters,options);

%% Evaluate best fit

probability = calculate_recapture_probability( ...
    releaseTimes, ...
    temperature, ...
    positionSamples, ...
    velocitySamples, ...
    physicalModel, ...
    options);

fittedCounts = background + amplitude*probability;

residuals = rawCounts-fittedCounts;
normalizedResiduals = residuals./rawCountErrors;

degreesOfFreedom = max(numel(rawCounts)-3,1);

chiSquare = sum(normalizedResiduals.^2);
reducedChiSquare = chiSquare/degreesOfFreedom;

totalSumSquares = sum((rawCounts-mean(rawCounts)).^2);

if totalSumSquares > 0
    rSquared = 1-sum(residuals.^2)/totalSumSquares;
else
    rSquared = NaN;
end

%% Smooth curve for plotting

plotTimes = linspace( ...
    min(releaseTimes), ...
    max(releaseTimes), ...
    options.NumFitCurvePoints).';

plotProbability = calculate_recapture_probability( ...
    plotTimes, ...
    temperature, ...
    positionSamples, ...
    velocitySamples, ...
    physicalModel, ...
    options);

plotCounts = background + amplitude*plotProbability;

%% Store result

fitResult.temperatureK = temperature;
fitResult.temperatureMicroK = temperature*1e6;

fitResult.amplitude = amplitude;
fitResult.background = background;

fitResult.releaseTimesSeconds = releaseTimes;
fitResult.rawCounts = rawCounts;
fitResult.rawCountErrors = rawCountErrors;

fitResult.recaptureProbability = probability;
fitResult.fittedCounts = fittedCounts;

fitResult.plotTimesSeconds = plotTimes;
fitResult.plotProbability = plotProbability;
fitResult.plotCounts = plotCounts;

fitResult.residuals = residuals;
fitResult.normalizedResiduals = normalizedResiduals;

fitResult.chiSquare = chiSquare;
fitResult.reducedChiSquare = reducedChiSquare;
fitResult.degreesOfFreedom = degreesOfFreedom;
fitResult.rSquared = rSquared;

fitResult.objectiveValue = objectiveValue;
fitResult.exitFlag = exitFlag;
fitResult.optimizerOutput = optimizerOutput;

end

%% Minimize Residuals
% Simulate the Montecarlo for given random positions and velocities
% And output the parameter sequence with smallest residuals
function objectiveValue = raw_count_objective( ...
parameters, ...
releaseTimes, ...
rawCounts, ...
rawCountErrors, ...
positionSamples, ...
velocitySamples, ...
physicalModel, ...
options)

[temperature,amplitude,background] = ...
    decode_fit_parameters(parameters,options);

probability = calculate_recapture_probability( ...
    releaseTimes, ...
    temperature, ...
    positionSamples, ...
    velocitySamples, ...
    physicalModel, ...
    options);

modelCounts = background + amplitude*probability;

residuals = ...
    (rawCounts-modelCounts)./rawCountErrors;

objectiveValue = sum(residuals.^2);

end

%% Function to help read out the parameters
function [temperature,amplitude,background] = ...
decode_fit_parameters(parameters,options)

% Logistic mapping keeps temperature inside the requested range.
fraction = 1/(1+exp(-parameters(1)));

temperature = ...
    options.MinimumTemperature + ...
    (options.MaximumTemperature-options.MinimumTemperature)*fraction;

% Exponential mapping keeps amplitude positive.
amplitude = exp(parameters(2));

background = parameters(3);

end

%% ========================================================================
% MONTE CARLO RECAPTURE MODEL
% ========================================================================

function probability = calculate_recapture_probability( ...
releaseTimes, ...
temperature, ...
positionSamples, ...
velocitySamples, ...
physicalModel, ...
options)

kB = physicalModel.boltzmannConstant;
mass = physicalModel.mass;

% Thermal position widths in the harmonic trap.
sigmaPosition = [ ...
    sqrt(kB*temperature/(mass*physicalModel.omegaX^2))
    sqrt(kB*temperature/(mass*physicalModel.omegaY^2))
    sqrt(kB*temperature/(mass*physicalModel.omegaZ^2)) ];

% Thermal velocity width.
sigmaVelocity = sqrt(kB*temperature/mass);

initialPositions = ...
    positionSamples .* sigmaPosition';

initialVelocities = ...
    velocitySamples*sigmaVelocity;

probability = zeros(numel(releaseTimes),1);

for i = 1:numel(releaseTimes)

    t = releaseTimes(i);

    % Ballistic expansion during release.
    positions = initialPositions + initialVelocities*t;
    velocities = initialVelocities;

    % Add gravity.
    axis = lower(options.GravityAxis);

    if axis == 'x'
        positions(:,1) = positions(:,1) - 0.5*options.Gravity*t^2;
        velocities(:,1) = velocities(:,1) - options.Gravity*t;

    elseif axis == 'y'
        positions(:,2) = positions(:,2) - 0.5*options.Gravity*t^2;
        velocities(:,2) = velocities(:,2) - options.Gravity*t;

    elseif axis == 'z'
        positions(:,3) = positions(:,3) - 0.5*options.Gravity*t^2;
        velocities(:,3) = velocities(:,3) - options.Gravity*t;
    end

    % Kinetic energy at the end of the release.
    kineticEnergy = ...
        0.5*mass*sum(velocities.^2,2);

    % Position-dependent trap depth.
    depthX = physicalModel.trapDepthX * ...
        exp(-2*positions(:,1).^2/options.WaistX^2);

    depthY = physicalModel.trapDepthY * ...
        exp(-2*positions(:,2).^2/options.WaistY^2);

    depthZ = physicalModel.trapDepthZ * ...
        exp(-2*positions(:,3).^2/options.WaistZ^2);

    trapDepth = depthX + depthY + depthZ;

    % A particle is recaptured if its kinetic energy is below
    % the local trap depth.
    probability(i) = mean(kineticEnergy < trapDepth);
end

end

%% ========================================================================
% COMBINE TWEEZER DATA
% ========================================================================

function [combinedX,combinedY,combinedError] = ...
combine_tweezer_scans(xCells,yCells,errorCells)

% Gather every release time used by any tweezer.
allX = vertcat(xCells{:});

if isempty(allX)
    combinedX = [];
    combinedY = [];
    combinedError = [];
    return;
end

combinedX = unique(allX);

combinedY = nan(size(combinedX));
combinedError = nan(size(combinedX));

for i = 1:numel(combinedX)

    x = combinedX(i);

    values = [];
    errors = [];

    for tweezer = 1:numel(xCells)

        if isempty(xCells{tweezer})
            continue;
        end

        % The release-time grids should normally match exactly.
        % The small tolerance makes this robust to tiny floating-point
        % differences.
        tolerance = max(1e-12,1e-9*max(1,abs(x)));

        index = find( ...
            abs(xCells{tweezer}-x) <= tolerance, ...
            1);

        if isempty(index)
            continue;
        end

        values(end+1,1) = yCells{tweezer}(index); %#ok<AGROW>
        errors(end+1,1) = errorCells{tweezer}(index); %#ok<AGROW>
    end

    [combinedY(i),combinedError(i)] = ...
        combine_mean_values(values,errors);
end

valid = ...
    isfinite(combinedX) & ...
    isfinite(combinedY) & ...
    isfinite(combinedError) & ...
    combinedError > 0;

combinedX = combinedX(valid);
combinedY = combinedY(valid);
combinedError = combinedError(valid);

end

function [meanValue,totalError] = combine_mean_values(values,errors)

values = values(:);
errors = errors(:);

valid = isfinite(values);

values = values(valid);
errors = errors(valid);

if isempty(values)
    meanValue = NaN;
    totalError = NaN;
    return;
end

meanValue = mean(values);
n = numel(values);

% Error on the mean from the individual tweezer uncertainties.
propagatedError = sqrt(sum(errors.^2))/n;

% Error on the mean from tweezer-to-tweezer scatter.
if n > 1
    scatterError = std(values)/sqrt(n);
else
    scatterError = 0;
end

% Treat the two sources of uncertainty as independent.
totalError = sqrt(propagatedError^2 + scatterError^2);

% This only matters for the special case of perfectly identical
% values with zero uncertainty.
if totalError == 0
    totalError = eps(max(abs(meanValue),1));
end

end

%% ========================================================================
% PLOTTING
% ========================================================================

function plot_single_fit( ...
xData,yData,yError,fitResult,scanID,tweezerNum,options,isCombined)

figure;
hold on;

dataHandle = errorbar( ...
    xData,yData,yError, ...
    'o','LineStyle','none', ...
    'DisplayName','Averaged raw counts');

fitX = fitResult.plotTimesSeconds / options.TimeScaleToSeconds;

fitHandle = plot( ...
    fitX,fitResult.plotCounts, ...
    '-','LineWidth',1.5, ...
    'DisplayName','Monte Carlo fit');

backgroundHandle = yline( ...
    fitResult.background,'--', ...
    'DisplayName','Fitted background');

xlabel(options.XAxisLabel,'Interpreter','none');
ylabel(options.YAxisLabel,'Interpreter','none');

set(gca,'XScale',options.XAxisScale);
set(gca,'YScale',options.YAxisScale);

if isCombined
    title(sprintf( ...
        '%s: all-tweezer average, scan ID %g', ...
        options.FitTitle,scanID), ...
        'Interpreter','none');
else
    title(sprintf( ...
        '%s: tweezer ROI %d, scan ID %g', ...
        options.FitTitle,tweezerNum,scanID), ...
        'Interpreter','none');
end

legend( ...
    [dataHandle fitHandle backgroundHandle], ...
    'Location','best');

grid on;
box on;

if options.ShowFitAnnotation

    text = { ...
        sprintf('T = %.4g microK',fitResult.temperatureMicroK)
        sprintf('Amplitude = %.5g',fitResult.amplitude)
        sprintf('Background = %.5g',fitResult.background)
        sprintf('\\chi^2_\\nu = %.4g',fitResult.reducedChiSquare)
        sprintf('R^2 = %.4g',fitResult.rSquared) };

    annotation( ...
        'textbox',[0.16 0.64 0.27 0.22], ...
        'String',text, ...
        'FitBoxToText','on', ...
        'BackgroundColor','white');
end

hold off;

end

function plot_all_combined_fits( ...
combinedX,combinedY,combinedYerr, ...
combinedFits,scanIDs,analyVar,options)

figure;
hold on;

handles = [];
labels = {};

for scan = 1:numel(scanIDs)

    if isempty(combinedFits{scan})
        continue;
    end

    color = analyVar.COLORS( ...
        mod(scan-1,size(analyVar.COLORS,1))+1,:);

    handles(end+1) = errorbar( ... %#ok<AGROW>
        combinedX{scan}, ...
        combinedY{scan}/max(combinedFits{scan}.plotCounts), ...
        combinedYerr{scan}/max(combinedFits{scan}.plotCounts), ...
        'o', ...
        'MarkerFaceColor',color, ...
        'MarkerEdgeColor','k', ...
        'Color',color);

    fitX = combinedFits{scan}.plotTimesSeconds / ...
        options.TimeScaleToSeconds;

    fitHandle = plot( ...
        fitX, ...
        combinedFits{scan}.plotCounts/max(combinedFits{scan}.plotCounts), ...
        '-', ...
        'LineWidth',1.5, ...
        'Color',color);

    fitHandle.HandleVisibility = 'off';

    labels{end+1} = sprintf('Scan ID %g',scanIDs(scan)); %#ok<AGROW>
end

xlabel(options.XAxisLabel,'Interpreter','none');
ylabel(options.YAxisLabel,'Interpreter','none');

title( ...
    'All-Tweezer Averaged Release-and-Recapture Fits', ...
    'Interpreter','none');

set(gca,'XScale',options.XAxisScale);
set(gca,'YScale',options.YAxisScale);

legend(handles,labels,'Location','best');

grid on;
box on;
hold off;

end



