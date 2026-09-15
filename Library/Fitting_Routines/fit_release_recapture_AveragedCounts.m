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
options.OmegaXHz = 50e3;
options.OmegaYHz = 50e3;
options.OmegaZHz = 10e3;

options.TrapWavelength = 532e-9;

options.WaistX = 1000e-9;
options.WaistY = 1000e-9;
options.WaistZ = 1000e-9;

options.Gravity = 9.8;
options.GravityAxis = 'z';

% Monte Carlo
options.NumParticles = 1000;
options.RandomSeed = 1;

% Fit parameters
options.InitialTemperature = 2e-6;
options.MinimumTemperature = 0.1e-6;
options.MaximumTemperature = 1000e-6;

options.InitialAmplitude = [];
options.InitialBackground = [];

options.MinimumError = 1;
options.NumFitCurvePoints = 50000;

options.FitOptions = optimset( ...
    'Display','off', ...
    'MaxFunEvals',3000, ...
    'MaxIter',1000, ...
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



% function funcOut = fit_release_recapture_AveragedCounts( ...
%     analyVar, indivDataset, avgDataset, useroptions)
% %FIT_RELEASE_RECAPTURE_AVERAGEDCOUNTS
% % Fit release-and-recapture measurements using averaged raw counts.
% %
% % The function:
% %   1. Uses get_averages for each tweezer ROI.
% %   2. Fits every averaged scan ID for every tweezer.
% %   3. Plots the averaged data and Monte Carlo fit for each tweezer.
% %   4. Averages corresponding tweezer points using both the uncertainty of
% %      each tweezer mean and the observed tweezer-to-tweezer spread.
% %   5. Fits and plots the properly combined all-tweezer average.
% %
% % Raw-count model:
% %
% %   Counts(t) = Background + Amplitude*P_recapture(t,T)
% %
% % INPUTS
% %   analyVar, indivDataset, avgDataset
% %       Standard analysis structures used by get_averages.
% %
% %   useroptions
% %       Optional structure. Supported fields include:
% %
% %       Analysis fields
% %       ---------------
% %       .IndVarField            Default: 'imagevcoAtom'
% %       .DepVarField            Default: 'OD_TotalCounts'
% %       .Statistics             Default: 'gaussian'
% %       .XAxisLabel             Default: IndVarField
% %       .YAxisLabel             Default: DepVarField
% %       .XAxisScale             Default: 'linear'
% %       .YAxisScale             Default: 'linear'
% %       .TimeScaleToSeconds     Multiply x data by this to obtain seconds.
% %                               Default: 0.001
% %
% %       Plot controls
% %       -------------
% %       .PlotIndividualTweezers Default: true
% %       .PlotAllTweezerAverage  Default: true
% %       .PlotAllFitsTogether    Default: false
% %       .ShowFitAnnotation      Default: true
% %       .FitTitle               Default:
% %                               'Release-and-Recapture Raw-Count Fit'
% %
% %       Monte Carlo parameters
% %       ----------------------
% %       .MassAMU                Default: 87.9056 Sr 88
% %       .OmegaXHz               Default: 50e3
% %       .OmegaYHz               Default: 50e3
% %       .OmegaZHz               Default: 10e3
% %       .TrapWavelength         Default: 532e-9 m
% %       .WaistX                 Default: 1000e-9 m
% %       .WaistY                 Default: 1000e-9 m
% %       .WaistZ                 Default: 1000e-9 m
% %       .Gravity                Default: 9.8 m/s^2
% %       .GravityAxis            Default: 'z'
% %       .NumParticles           Default: 100000
% %       .RandomSeed             Default: 1
% %
% %       Fit controls
% %       ------------
% %       .InitialTemperature     Default: 30e-6 K
% %       .MinimumTemperature     Default: 0.1e-6 K
% %       .MaximumTemperature     Default: 1000e-6 K
% %       .InitialAmplitude       Default: determined from each data set
% %       .InitialBackground      Default: determined from each data set
% %       .MinimumError           Default: 1
% %       .NumFitCurvePoints      Default: 500
% %       .FitOptions             Default: optimset('Display','off')
% %
% % OUTPUT
% %   results
% %       Structure containing individual-tweezer and combined-tweezer fits.
% 
%     if nargin < 4 || isempty(useroptions)
%         useroptions = struct;
%     end
% 
%     options = create_default_options();
% 
%     optionNames = fieldnames(useroptions);
% 
%     for optionIndex = 1:numel(optionNames)
%         optionName = optionNames{optionIndex};
%         options.(optionName) = useroptions.(optionName);
%     end
% 
%     if isempty(options.XAxisLabel)
%         options.XAxisLabel = options.IndVarField;
%     end
% 
%     if isempty(options.YAxisLabel)
%         options.YAxisLabel = options.DepVarField;
%     end
% 
%     scanIDs = analyVar.uniqScanList(:);
%     numScanIDs = numel(scanIDs);
% 
%     if analyVar.UseTweezer
%         roiFile = fullfile(analyVar.analyOutDir,'tweezerROI.mat');
%         roiData = load(roiFile,'tweezerROI');
%         numTweezers = size(roiData.tweezerROI.centersXY,1);
%     else
%         numTweezers = 1;
%     end
% 
%     fprintf('Fitting release-and-recapture raw-count data.\n');
%     fprintf('Number of tweezer ROIs: %d\n',numTweezers);
%     fprintf('Number of scan IDs: %d\n\n',numScanIDs);
% 
%     %% Fixed Monte Carlo samples
% 
%     previousRandomState = rng;
%     randomStateCleanup = onCleanup(@() rng(previousRandomState)); %#ok<NASGU>
% 
%     rng(options.RandomSeed,'twister');
% 
%     standardPositionSamples = randn(options.NumParticles,3);
%     standardVelocitySamples = randn(options.NumParticles,3);
% 
%     physicalModel = create_physical_model(options);
% 
%     %% Storage for each tweezer
% 
%     xByTweezer = cell(numTweezers,numScanIDs);
%     yByTweezer = cell(numTweezers,numScanIDs);
%     yerrByTweezer = cell(numTweezers,numScanIDs);
% 
%     individualFits = cell(numTweezers,numScanIDs);
% 
%     %% Average and fit each tweezer separately
% 
%     for tweezerNum = 1:numTweezers
% 
%         if analyVar.UseTweezer
%             fprintf('Processing tweezer ROI %d of %d.\n', ...
%                 tweezerNum,numTweezers);
%         end
% 
%         [xavg,yavg,yerr] = get_averages( ...
%             analyVar, ...
%             indivDataset, ...
%             avgDataset, ...
%             options.IndVarField, ...
%             options.DepVarField, ...
%             options.Statistics, ...
%             tweezerNum);
% 
%         for scanIndex = 1:numScanIDs
% 
%             [xData,yData,errorData] = clean_scan_data( ...
%                 xavg{scanIndex}, ...
%                 yavg{scanIndex}, ...
%                 yerr{scanIndex}, ...
%                 options.MinimumError);
% 
%             xByTweezer{tweezerNum,scanIndex} = xData;
%             yByTweezer{tweezerNum,scanIndex} = yData;
%             yerrByTweezer{tweezerNum,scanIndex} = errorData;
% 
%             if numel(xData) < 3
%                 warning(['Tweezer %d, scan ID %g has fewer than three ' ...
%                          'valid points and will not be fitted.'], ...
%                     tweezerNum,scanIDs(scanIndex));
%                 continue;
%             end
% 
%             releaseTimesSeconds = ...
%                 xData*options.TimeScaleToSeconds;
% 
%             thisFit = fit_single_raw_count_scan( ...
%                 releaseTimesSeconds, ...
%                 yData, ...
%                 errorData, ...
%                 standardPositionSamples, ...
%                 standardVelocitySamples, ...
%                 physicalModel, ...
%                 options);
% 
%             thisFit.xDataOriginalUnits = xData;
%             thisFit.scanID = scanIDs(scanIndex);
%             thisFit.tweezerNum = tweezerNum;
% 
%             individualFits{tweezerNum,scanIndex} = thisFit;
% 
%             if options.PlotIndividualTweezers
%                 plot_single_fit( ...
%                     xData, ...
%                     yData, ...
%                     errorData, ...
%                     thisFit, ...
%                     scanIDs(scanIndex), ...
%                     tweezerNum, ...
%                     options, ...
%                     false);
%             end
%         end
%     end
% 
%     %% Properly combine corresponding tweezer averages
% 
%     combinedX = cell(numScanIDs,1);
%     combinedY = cell(numScanIDs,1);
%     combinedYerr = cell(numScanIDs,1);
%     combinedFits = cell(numScanIDs,1);
% 
%     for scanIndex = 1:numScanIDs
% 
%         [combinedX{scanIndex}, ...
%          combinedY{scanIndex}, ...
%          combinedYerr{scanIndex}] = combine_tweezer_scans( ...
%             xByTweezer(:,scanIndex), ...
%             yByTweezer(:,scanIndex), ...
%             yerrByTweezer(:,scanIndex));
% 
%         if numel(combinedX{scanIndex}) < 3
%             warning(['Combined tweezer data for scan ID %g has fewer ' ...
%                      'than three valid points and will not be fitted.'], ...
%                 scanIDs(scanIndex));
%             continue;
%         end
% 
%         releaseTimesSeconds = ...
%             combinedX{scanIndex}*options.TimeScaleToSeconds;
% 
%         thisFit = fit_single_raw_count_scan( ...
%             releaseTimesSeconds, ...
%             combinedY{scanIndex}, ...
%             combinedYerr{scanIndex}, ...
%             standardPositionSamples, ...
%             standardVelocitySamples, ...
%             physicalModel, ...
%             options);
% 
%         thisFit.xDataOriginalUnits = combinedX{scanIndex};
%         thisFit.scanID = scanIDs(scanIndex);
%         thisFit.tweezerNum = NaN;
% 
%         combinedFits{scanIndex} = thisFit;
% 
%         if options.PlotAllTweezerAverage
%             plot_single_fit( ...
%                 combinedX{scanIndex}, ...
%                 combinedY{scanIndex}, ...
%                 combinedYerr{scanIndex}, ...
%                 thisFit, ...
%                 scanIDs(scanIndex), ...
%                 NaN, ...
%                 options, ...
%                 true);
%         end
%     end
% 
%     %% Optional overview plot of all combined scans and fits
% 
%     if options.PlotAllFitsTogether
%         plot_all_combined_fits( ...
%             combinedX, ...
%             combinedY, ...
%             combinedYerr, ...
%             combinedFits, ...
%             scanIDs, ...
%             analyVar, ...
%             options);
%     end
% 
%     %% Return results
% 
%     %   results = struct;
% %   
%     %   results.options = options;
%     %   results.scanIDs = scanIDs;
%     %   results.numTweezers = numTweezers;
% %   
%     %   results.individual.x = xByTweezer;
%     %   results.individual.y = yByTweezer;
%     %   results.individual.yerr = yerrByTweezer;
%     %   results.individual.fit = individualFits;
% %   
%     %   results.combined.x = combinedX;
%     %   results.combined.y = combinedY;
%     %   results.combined.yerr = combinedYerr;
%     %   results.combined.fit = combinedFits;
% %   
%     %   results.physicalModel = physicalModel;
% 
%     funcOut.analyVar = analyVar;
%     funcOut.indivDataset = indivDataset;
%     funcOut.avgDataset = avgDataset;
% 
% end
% 
% 
% %% ======================================================================
% %                   OPTIONS
% % =======================================================================
% 
% function options = create_default_options()
% 
%     options = struct;
% 
%     options.IndVarField = 'imagevcoAtom';
%     options.DepVarField = 'OD_TotalCounts';
%     options.Statistics = 'gaussian';
% 
%     options.XAxisLabel = '';
%     options.YAxisLabel = '';
% 
%     options.XAxisScale = 'linear';
%     options.YAxisScale = 'linear';
% 
%     options.TimeScaleToSeconds = 0.001; %% Labview records in ms
% 
%     options.PlotIndividualTweezers = false;
%     options.PlotAllTweezerAverage = true;
%     options.PlotAllFitsTogether = true;
%     options.ShowFitAnnotation = true;
% 
%     options.FitTitle = ...
%         'Release-and-Recapture Raw-Count Fit';
% 
%     options.MassAMU = 87.9056; %% Sr88
% 
%     options.OmegaXHz = 50e3;
%     options.OmegaYHz = 50e3;
%     options.OmegaZHz = 10e3;
% 
%     options.TrapWavelength = 532e-9;
% 
%     options.WaistX = 1000e-9;
%     options.WaistY = 1000e-9;
%     options.WaistZ = 1000e-9;
% 
%     options.Gravity = 9.8;
%     options.GravityAxis = 'z';
% 
%     options.NumParticles = 1000;
%     options.RandomSeed = 1;
% 
%     options.InitialTemperature = 10e-6;
%     options.MinimumTemperature = 0.1e-6;
%     options.MaximumTemperature = 1000e-6;
% 
%     options.InitialAmplitude = [];
%     options.InitialBackground = [];
% 
%     options.MinimumError = 1;
%     options.NumFitCurvePoints = 500;
% 
%     options.FitOptions = optimset( ...
%         'Display','off', ...
%         'MaxFunEvals',3000, ...
%         'MaxIter',1000, ...
%         'TolX',1e-7, ...
%         'TolFun',1e-7);
% 
% end
% 
% 
% function physicalModel = create_physical_model(options)
% 
%     boltzmannConstant = 1.380649e-23;
%     atomicMassUnit = 1.66053906660e-27;
% 
%     mass = options.MassAMU*atomicMassUnit;
% 
%     omegaX = 2*pi*options.OmegaXHz;
%     omegaY = 2*pi*options.OmegaYHz;
%     omegaZ = 2*pi*options.OmegaZHz;
% 
%     waveNumber = 2*pi/options.TrapWavelength;
% 
%     physicalModel = struct;
% 
%     physicalModel.boltzmannConstant = boltzmannConstant;
%     physicalModel.mass = mass;
% 
%     physicalModel.omegaX = omegaX;
%     physicalModel.omegaY = omegaY;
%     physicalModel.omegaZ = omegaZ;
% 
%     physicalModel.trapDepthX = ...
%         mass*omegaX^2/(2*waveNumber^2);
% 
%     physicalModel.trapDepthY = ...
%         mass*omegaY^2/(2*waveNumber^2);
% 
%     physicalModel.trapDepthZ = ...
%         mass*omegaZ^2/(2*waveNumber^2);
% 
% end
% 
% 
% function [xData,yData,errorData] = ...
%     clean_scan_data(xData,yData,errorData,minimumError)
% 
%     xData = xData(:);
%     yData = yData(:);
% 
%     if isempty(errorData)
%         errorData = nan(size(yData));
%     else
%         errorData = errorData(:);
%     end
% 
%     if numel(xData) ~= numel(yData)
%         error('The x and y data must have the same number of elements.');
%     end
% 
%     if numel(errorData) ~= numel(yData)
%         error('The y data and y-error data must have the same length.');
%     end
% 
%     valid = isfinite(xData) & isfinite(yData);
% 
%     xData = xData(valid);
%     yData = yData(valid);
%     errorData = errorData(valid);
% 
%     [xData,sortOrder] = sort(xData);
%     yData = yData(sortOrder);
%     errorData = errorData(sortOrder);
% 
%     validPositiveErrors = ...
%         errorData(isfinite(errorData) & errorData > 0);
% 
%     if isempty(validPositiveErrors)
%         replacementError = minimumError;
%     else
%         replacementError = max( ...
%             median(validPositiveErrors), ...
%             minimumError);
%     end
% 
%     invalidErrors = ...
%         ~isfinite(errorData) | errorData <= 0;
% 
%     errorData(invalidErrors) = replacementError;
%     errorData = max(errorData,minimumError);
% 
% end
% 
% 
% function fitResult = fit_single_raw_count_scan( ...
%     releaseTimes, ...
%     rawCounts, ...
%     rawCountErrors, ...
%     standardPositionSamples, ...
%     standardVelocitySamples, ...
%     physicalModel, ...
%     options)
% 
%     if isempty(options.InitialAmplitude)
%         initialAmplitude = max(rawCounts)-min(rawCounts);
% 
%         if ~isfinite(initialAmplitude) || initialAmplitude <= 0
%             initialAmplitude = max(abs(rawCounts));
% 
%             if ~isfinite(initialAmplitude) || initialAmplitude <= 0
%                 initialAmplitude = 1;
%             end
%         end
%     else
%         initialAmplitude = options.InitialAmplitude;
%     end
% 
%     if isempty(options.InitialBackground)
%         initialBackground = min(rawCounts);
%     else
%         initialBackground = options.InitialBackground;
%     end
% 
%     temperatureFraction = ...
%         (options.InitialTemperature-options.MinimumTemperature) / ...
%         (options.MaximumTemperature-options.MinimumTemperature);
% 
%     temperatureFraction = min( ...
%         max(temperatureFraction,1e-9), ...
%         1-1e-9);
% 
%     initialTemperatureParameter = ...
%         log(temperatureFraction/(1-temperatureFraction));
% 
%     initialParameters = [
%         initialTemperatureParameter
%         log(max(initialAmplitude,eps))
%         initialBackground
%     ];
% 
%     objectiveFunction = @(fitParameters) ...
%         raw_count_objective( ...
%             fitParameters, ...
%             releaseTimes, ...
%             rawCounts, ...
%             rawCountErrors, ...
%             standardPositionSamples, ...
%             standardVelocitySamples, ...
%             physicalModel, ...
%             options);
% 
%     [bestParameters,objectiveValue,exitFlag,optimizerOutput] = ...
%         fminsearch( ...
%             objectiveFunction, ...
%             initialParameters, ...
%             options.FitOptions);
% 
%     [temperature,amplitude,background] = ...
%         decode_fit_parameters(bestParameters,options);
% 
%     recaptureProbability = calculate_recapture_probability( ...
%         releaseTimes, ...
%         temperature, ...
%         standardPositionSamples, ...
%         standardVelocitySamples, ...
%         physicalModel, ...
%         options);
% 
%     fittedCounts = ...
%         background + amplitude*recaptureProbability;
% 
%     residuals = rawCounts-fittedCounts;
%     normalizedResiduals = residuals./rawCountErrors;
% 
%     numFitParameters = 3;
%     degreesOfFreedom = max(numel(rawCounts)-numFitParameters,1);
% 
%     chiSquare = sum(normalizedResiduals.^2);
%     reducedChiSquare = chiSquare/degreesOfFreedom;
% 
%     totalSumSquares = ...
%         sum((rawCounts-mean(rawCounts)).^2);
% 
%     if totalSumSquares > 0
%         rSquared = ...
%             1-sum(residuals.^2)/totalSumSquares;
%     else
%         rSquared = NaN;
%     end
% 
%     plotTimes = linspace( ...
%         min(releaseTimes), ...
%         max(releaseTimes), ...
%         options.NumFitCurvePoints).';
% 
%     plotProbability = calculate_recapture_probability( ...
%         plotTimes, ...
%         temperature, ...
%         standardPositionSamples, ...
%         standardVelocitySamples, ...
%         physicalModel, ...
%         options);
% 
%     plotCounts = ...
%         background + amplitude*plotProbability;
% 
%     fitResult = struct;
% 
%     fitResult.temperatureK = temperature;
%     fitResult.temperatureMicroK = temperature*1e6;
% 
%     fitResult.amplitude = amplitude;
%     fitResult.background = background;
% 
%     fitResult.releaseTimesSeconds = releaseTimes;
%     fitResult.rawCounts = rawCounts;
%     fitResult.rawCountErrors = rawCountErrors;
% 
%     fitResult.recaptureProbability = recaptureProbability;
%     fitResult.fittedCounts = fittedCounts;
% 
%     fitResult.plotTimesSeconds = plotTimes;
%     fitResult.plotProbability = plotProbability;
%     fitResult.plotCounts = plotCounts;
% 
%     fitResult.residuals = residuals;
%     fitResult.normalizedResiduals = normalizedResiduals;
% 
%     fitResult.chiSquare = chiSquare;
%     fitResult.reducedChiSquare = reducedChiSquare;
%     fitResult.degreesOfFreedom = degreesOfFreedom;
%     fitResult.rSquared = rSquared;
% 
%     fitResult.objectiveValue = objectiveValue;
%     fitResult.exitFlag = exitFlag;
%     fitResult.optimizerOutput = optimizerOutput;
% 
% end
% 
% 
% function objectiveValue = raw_count_objective( ...
%     fitParameters, ...
%     releaseTimes, ...
%     rawCounts, ...
%     rawCountErrors, ...
%     standardPositionSamples, ...
%     standardVelocitySamples, ...
%     physicalModel, ...
%     options)
% 
%     [temperature,amplitude,background] = ...
%         decode_fit_parameters(fitParameters,options);
% 
%     probability = calculate_recapture_probability( ...
%         releaseTimes, ...
%         temperature, ...
%         standardPositionSamples, ...
%         standardVelocitySamples, ...
%         physicalModel, ...
%         options);
% 
%     modelCounts = ...
%         background + amplitude*probability;
% 
%     normalizedResiduals = ...
%         (rawCounts-modelCounts)./rawCountErrors;
% 
%     objectiveValue = sum(normalizedResiduals.^2);
% 
%     if ~isfinite(objectiveValue)
%         objectiveValue = realmax;
%     end
% 
% end
% 
% 
% function [temperature,amplitude,background] = ...
%     decode_fit_parameters(fitParameters,options)
% 
%     temperatureFraction = ...
%         1/(1+exp(-fitParameters(1)));
% 
%     temperature = ...
%         options.MinimumTemperature + ...
%         (options.MaximumTemperature-options.MinimumTemperature) * ...
%         temperatureFraction;
% 
%     amplitude = exp(fitParameters(2));
%     background = fitParameters(3);
% 
% end
% 
% 
% function recaptureProbability = calculate_recapture_probability( ...
%     releaseTimes, ...
%     temperature, ...
%     standardPositionSamples, ...
%     standardVelocitySamples, ...
%     physicalModel, ...
%     options)
% 
%     kb = physicalModel.boltzmannConstant;
%     mass = physicalModel.mass;
% 
%     positionStandardDeviations = [
%         sqrt(kb*temperature/(mass*physicalModel.omegaX^2))
%         sqrt(kb*temperature/(mass*physicalModel.omegaY^2))
%         sqrt(kb*temperature/(mass*physicalModel.omegaZ^2))
%     ];
% 
%     velocityStandardDeviation = ...
%         sqrt(kb*temperature/mass);
% 
%     initialPositions = ...
%         standardPositionSamples .* ...
%         reshape(positionStandardDeviations,1,3);
% 
%     initialVelocities = ...
%         standardVelocitySamples*velocityStandardDeviation;
% 
%     numTimes = numel(releaseTimes);
%     numParticles = size(initialPositions,1);
% 
%     recaptureProbability = nan(numTimes,1);
% 
%     for timeIndex = 1:numTimes
% 
%         releaseTime = releaseTimes(timeIndex);
% 
%         positions = ...
%             initialPositions + initialVelocities*releaseTime;
% 
%         velocities = initialVelocities;
% 
%         switch lower(options.GravityAxis)
% 
%             case 'x'
%                 positions(:,1) = positions(:,1) ...
%                     - 0.5*options.Gravity*releaseTime^2;
% 
%                 velocities(:,1) = velocities(:,1) ...
%                     - options.Gravity*releaseTime;
% 
%             case 'y'
%                 positions(:,2) = positions(:,2) ...
%                     - 0.5*options.Gravity*releaseTime^2;
% 
%                 velocities(:,2) = velocities(:,2) ...
%                     - options.Gravity*releaseTime;
% 
%             case 'z'
%                 positions(:,3) = positions(:,3) ...
%                     - 0.5*options.Gravity*releaseTime^2;
% 
%                 velocities(:,3) = velocities(:,3) ...
%                     - options.Gravity*releaseTime;
%         end
% 
%         kineticEnergy = ...
%             0.5*mass*sum(velocities.^2,2);
% 
%         localTrapDepthX = ...
%             physicalModel.trapDepthX .* ...
%             exp(-2*positions(:,1).^2/options.WaistX^2);
% 
%         localTrapDepthY = ...
%             physicalModel.trapDepthY .* ...
%             exp(-2*positions(:,2).^2/options.WaistY^2);
% 
%         localTrapDepthZ = ...
%             physicalModel.trapDepthZ .* ...
%             exp(-2*positions(:,3).^2/options.WaistZ^2);
% 
%         localTrapDepth = ...
%             localTrapDepthX + ...
%             localTrapDepthY + ...
%             localTrapDepthZ;
% 
%         recaptured = kineticEnergy < localTrapDepth;
% 
%         recaptureProbability(timeIndex) = ...
%             sum(recaptured)/numParticles;
%     end
% 
% end
% 
% 
% function [combinedX,combinedY,combinedError] = ...
%     combine_tweezer_scans(xCells,yCells,errorCells)
% 
%     numTweezers = numel(xCells);
% 
%     allX = [];
% 
%     for tweezerNum = 1:numTweezers
%         if ~isempty(xCells{tweezerNum})
%             allX = [allX; xCells{tweezerNum}(:)]; %#ok<AGROW>
%         end
%     end
% 
%     if isempty(allX)
%         combinedX = [];
%         combinedY = [];
%         combinedError = [];
%         return;
%     end
% 
%     combinedX = unique(allX,'sorted');
% 
%     numX = numel(combinedX);
% 
%     combinedY = nan(numX,1);
%     combinedError = nan(numX,1);
% 
%     for xIndex = 1:numX
% 
%         currentX = combinedX(xIndex);
%         xTolerance = max(1e-12,1e-9*max(1,abs(currentX)));
% 
%         values = [];
%         errors = [];
% 
%         for tweezerNum = 1:numTweezers
% 
%             thisX = xCells{tweezerNum};
%             thisY = yCells{tweezerNum};
%             thisError = errorCells{tweezerNum};
% 
%             if isempty(thisX)
%                 continue;
%             end
% 
%             matchingIndex = find( ...
%                 abs(thisX-currentX) <= xTolerance, ...
%                 1, ...
%                 'first');
% 
%             if isempty(matchingIndex)
%                 continue;
%             end
% 
%             thisValue = thisY(matchingIndex);
%             thisUncertainty = thisError(matchingIndex);
% 
%             if isfinite(thisValue)
%                 values(end+1,1) = thisValue; %#ok<AGROW>
%                 errors(end+1,1) = thisUncertainty; %#ok<AGROW>
%             end
%         end
% 
%         [combinedY(xIndex),combinedError(xIndex)] = ...
%             combine_mean_values(values,errors);
%     end
% 
%     valid = ...
%         isfinite(combinedX) & ...
%         isfinite(combinedY) & ...
%         isfinite(combinedError) & ...
%         combinedError > 0;
% 
%     combinedX = combinedX(valid);
%     combinedY = combinedY(valid);
%     combinedError = combinedError(valid);
% 
% end
% 
% 
% function [combinedMean,combinedError] = ...
%     combine_mean_values(values,errors)
% 
%     values = values(:);
%     errors = errors(:);
% 
%     validValues = isfinite(values);
% 
%     values = values(validValues);
%     errors = errors(validValues);
% 
%     numValues = numel(values);
% 
%     if numValues == 0
%         combinedMean = NaN;
%         combinedError = NaN;
%         return;
%     end
% 
%     combinedMean = mean(values);
% 
%     validErrors = isfinite(errors) & errors >= 0;
% 
%     if any(validErrors)
%         propagatedError = ...
%             sqrt(sum(errors(validErrors).^2))/numValues;
%     else
%         propagatedError = 0;
%     end
% 
%     if numValues > 1
%         scatterError = ...
%             std(values,0)/sqrt(numValues);
%     else
%         scatterError = 0;
%     end
% 
%     combinedError = sqrt( ...
%         propagatedError^2 + scatterError^2);
% 
%     if combinedError <= 0 || ~isfinite(combinedError)
%         combinedError = eps(max(abs(combinedMean),1));
%     end
% 
% end
% 
% 
% function plot_single_fit( ...
%     xData, ...
%     yData, ...
%     yError, ...
%     fitResult, ...
%     scanID, ...
%     tweezerNum, ...
%     options, ...
%     isCombined)
% 
%     figure;
%     hold on;
% 
%     dataHandle = errorbar( ...
%         xData, ...
%         yData, ...
%         yError, ...
%         'o', ...
%         'LineStyle','none', ...
%         'DisplayName','Averaged raw counts');
% 
%     fitXOriginalUnits = ...
%         fitResult.plotTimesSeconds/options.TimeScaleToSeconds;
% 
%     fitHandle = plot( ...
%         fitXOriginalUnits, ...
%         fitResult.plotCounts, ...
%         '-', ...
%         'LineWidth',1.5, ...
%         'DisplayName','Monte Carlo fit');
% 
%     backgroundHandle = yline( ...
%         fitResult.background, ...
%         '--', ...
%         'DisplayName','Fitted background');
% 
%     xlabel(options.XAxisLabel,'Interpreter','none');
%     ylabel(options.YAxisLabel,'Interpreter','none');
% 
%     set(gca,'XScale',options.XAxisScale);
%     set(gca,'YScale',options.YAxisScale);
% 
%     if isCombined
%         titleText = sprintf( ...
%             '%s: all-tweezer average, scan ID %g', ...
%             options.FitTitle, ...
%             scanID);
%     else
%         titleText = sprintf( ...
%             '%s: tweezer ROI %d, scan ID %g', ...
%             options.FitTitle, ...
%             tweezerNum, ...
%             scanID);
%     end
% 
%     title(titleText,'Interpreter','none');
% 
%     legend( ...
%         [dataHandle fitHandle backgroundHandle], ...
%         'Location','best');
% 
%     grid on;
%     box on;
% 
%     if options.ShowFitAnnotation
% 
%         fitText = {
%             sprintf('T = %.4g microK', ...
%                 fitResult.temperatureMicroK)
%             sprintf('Amplitude = %.5g', ...
%                 fitResult.amplitude)
%             sprintf('Background = %.5g', ...
%                 fitResult.background)
%             sprintf('\\chi^2_\\nu = %.4g', ...
%                 fitResult.reducedChiSquare)
%             sprintf('R^2 = %.4g', ...
%                 fitResult.rSquared)
%         };
% 
%         annotation( ...
%             'textbox', ...
%             [0.16 0.64 0.27 0.22], ...
%             'String',fitText, ...
%             'FitBoxToText','on', ...
%             'BackgroundColor','white');
%     end
% 
%     hold off;
% 
% end
% 
% 
% function plot_all_combined_fits( ...
%     combinedX, ...
%     combinedY, ...
%     combinedYerr, ...
%     combinedFits, ...
%     scanIDs, ...
%     analyVar, ...
%     options)
% 
%     figure;
%     hold on;
% 
%     numScanIDs = numel(scanIDs);
%     legendHandles = gobjects(numScanIDs,1);
%     legendLabels = cell(numScanIDs,1);
%     validLegendEntry = false(numScanIDs,1);
% 
%     for scanIndex = 1:numScanIDs
% 
%         if isempty(combinedFits{scanIndex})
%             continue;
%         end
% 
%         colorIndex = mod(scanIndex-1,size(analyVar.COLORS,1))+1;
%         thisColor = analyVar.COLORS(colorIndex,:);
% 
%         legendHandles(scanIndex) = errorbar( ...
%             combinedX{scanIndex}, ...
%             combinedY{scanIndex}, ...
%             combinedYerr{scanIndex}, ...
%             'o', ...
%             'LineStyle','none', ...
%             'MarkerFaceColor',thisColor, ...
%             'MarkerEdgeColor','k', ...
%             'Color',thisColor);
% 
%         fitXOriginalUnits = ...
%             combinedFits{scanIndex}.plotTimesSeconds / ...
%             options.TimeScaleToSeconds;
% 
%         fitHandle = plot( ...
%             fitXOriginalUnits, ...
%             combinedFits{scanIndex}.plotCounts, ...
%             '-', ...
%             'LineWidth',1.5, ...
%             'Color',thisColor);
% 
%         fitHandle.HandleVisibility = 'off';
% 
%         legendLabels{scanIndex} = ...
%             sprintf('Scan ID %g',scanIDs(scanIndex));
% 
%         validLegendEntry(scanIndex) = true;
%     end
% 
%     xlabel(options.XAxisLabel,'Interpreter','none');
%     ylabel(options.YAxisLabel,'Interpreter','none');
% 
%     title( ...
%         'All-Tweezer Averaged Release-and-Recapture Fits', ...
%         'Interpreter','none');
% 
%     set(gca,'XScale',options.XAxisScale);
%     set(gca,'YScale',options.YAxisScale);
% 
%     legend( ...
%         legendHandles(validLegendEntry), ...
%         legendLabels(validLegendEntry), ...
%         'Location','best');
% 
%     grid on;
%     box on;
%     hold off;
% 
% end

