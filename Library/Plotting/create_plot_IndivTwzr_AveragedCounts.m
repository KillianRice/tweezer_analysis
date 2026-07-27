function create_plot_IndivTwzr_AveragedCounts( ...
    analyVar, indivDataset, avgDataset, plotRawCounts)

% Plot individually averaged ROI signals and averages over ROI groups.
%
% ROI ordering:
%   1 : numRealTweezers              Real tweezers
%   Following ROIs                   Non-tweezer section 1
%   Final ROIs                       Non-tweezer section 2
%
% For an unweighted mean of N ROI values, m_i, with uncertainties s_i:
%
%   meanValue = mean(m_i)
%
%   propagatedError = sqrt(sum(s_i.^2))/N
%
% The final error also includes the observed spread between ROI means:
%
%   scatterError = std(m_i)/sqrt(N)
%
%   totalError = sqrt(propagatedError^2 + scatterError^2)


%% Select variables

indVarField = 'imagevcoAtom';

if plotRawCounts
    depVarField = 'OD_TotalCountsImg1Raw';
else
    depVarField = 'OD_TotalCounts';
end

[xdataClean, ydataClean] = getxy( ...
    indVarField, ...
    depVarField, ...
    analyVar, ...
    indivDataset, ...
    avgDataset);

scanIDs = analyVar.uniqScanList(:);
numScanIDs = numel(scanIDs);


%% Determine number of ROIs

roiFile = fullfile( ...
    analyVar.analyOutDir, ...
    'tweezerROI.mat');

roiData = load(roiFile,'tweezerROI');
tweezerROI = roiData.tweezerROI;

numROIs = size(tweezerROI.centersXY,1);


%% Define ROI groups

numNonTweezers1 = analyVar.numFakeTweezers1;
numNonTweezers2 = analyVar.numFakeTweezers2;

validateattributes( ...
    numNonTweezers1, ...
    {'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});

validateattributes( ...
    numNonTweezers2, ...
    {'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});

numRealTweezers = ...
    numROIs - numNonTweezers1 - numNonTweezers2;

if numRealTweezers < 1
    error(['The requested non-tweezer groups leave no real ' ...
           'tweezer ROIs.']);
end

realTweezerIndices = ...
    1:numRealTweezers;

nonTweezer1Indices = ...
    (numRealTweezers + 1): ...
    (numRealTweezers + numNonTweezers1);

nonTweezer2Indices = ...
    (numRealTweezers + numNonTweezers1 + 1): ...
    numROIs;


%% Preallocate results

x = cell(numScanIDs,1);

y = cell(numScanIDs,numROIs);
yerr = cell(numScanIDs,numROIs);

% Scalar mean of each complete ROI curve
avgSignal = nan(numScanIDs,numROIs);
avgSignalErr = nan(numScanIDs,numROIs);

% Point-by-point group curves
yTweezerMean = cell(numScanIDs,1);
yTweezerErr = cell(numScanIDs,1);

yNonTweezer1Mean = cell(numScanIDs,1);
yNonTweezer1Err = cell(numScanIDs,1);

yNonTweezer2Mean = cell(numScanIDs,1);
yNonTweezer2Err = cell(numScanIDs,1);

% Scalar group averages
avgAllTweezers = nan(numScanIDs,1);
avgAllTweezersErr = nan(numScanIDs,1);

avgAllNonTweezers1 = nan(numScanIDs,1);
avgAllNonTweezers1Err = nan(numScanIDs,1);

avgAllNonTweezers2 = nan(numScanIDs,1);
avgAllNonTweezers2Err = nan(numScanIDs,1);

fprintf('Plotting individual tweezer counts.\n\n');


%% Process each scan ID

for id = 1:numScanIDs

    %% Collect all x values for this scan ID

    xAll = [];

    for basename = 1:analyVar.numBasenamesAtom

        if scanIDs(id) ~= analyVar.meanListVar(basename)
            continue;
        end

        thisX = xdataClean{basename}(:);
        xAll = [xAll; thisX]; %#ok<AGROW>
    end

    if isempty(xAll)
        warning( ...
            'No data found for scan ID %g.', ...
            scanIDs(id));

        continue;
    end

    x{id} = unique(xAll,'sorted');
    numX = numel(x{id});


    %% Initialize ROI curves for this scan ID

    for roiNum = 1:numROIs

        y{id,roiNum} = nan(numX,1);
        yerr{id,roiNum} = nan(numX,1);
    end


    %% Average repeated measurements at each x value

    for xIndex = 1:numX

        currentX = x{id}(xIndex);

        valuesAtX = cell(1,numROIs);

        for basename = 1:analyVar.numBasenamesAtom

            if scanIDs(id) ~= analyVar.meanListVar(basename)
                continue;
            end

            xBase = xdataClean{basename}(:);
            yBase = ydataClean{basename};

            if size(yBase,1) ~= numel(xBase)
                error(['The number of rows in ydataClean{%d} does not ' ...
                       'match the number of corresponding x values.'], ...
                    basename);
            end

            xTolerance = max( ...
                1e-12, ...
                1e-9*max(1,abs(currentX)));

            matchingRows = find( ...
                abs(xBase-currentX) <= xTolerance);

            numAvailableROIs = min( ...
                numROIs, ...
                size(yBase,2));

            for rowNum = matchingRows(:)'

                for roiNum = 1:numAvailableROIs

                    thisValue = yBase(rowNum,roiNum);

                    if isfinite(thisValue)

                        valuesAtX{roiNum}(end+1,1) = ...
                            thisValue;
                    end
                end
            end
        end


        %% Mean and standard error for each individual ROI

        for roiNum = 1:numROIs

            values = valuesAtX{roiNum};
            values = values(isfinite(values));

            numValues = numel(values);

            if numValues == 0
                continue;
            end

            y{id,roiNum}(xIndex) = mean(values);

            if numValues > 1

                yerr{id,roiNum}(xIndex) = ...
                    std(values,0) / sqrt(numValues);

            else

                yerr{id,roiNum}(xIndex) = 0;
            end
        end
    end


    %% Scalar average of the complete curve for each ROI

    for roiNum = 1:numROIs

        thisCurve = y{id,roiNum};
        thisCurveErr = yerr{id,roiNum};

        [avgSignal(id,roiNum), ...
         avgSignalErr(id,roiNum)] = ...
            combineMeanValues( ...
                thisCurve, ...
                thisCurveErr);
    end


    %% Point-by-point average of real-tweezer curves

    realCurveMatrix = ...
        cell2mat(y(id,realTweezerIndices));

    realErrorMatrix = ...
        cell2mat(yerr(id,realTweezerIndices));

    [yTweezerMean{id}, ...
     yTweezerErr{id}] = ...
        combineMeanMatrix( ...
            realCurveMatrix, ...
            realErrorMatrix);


    %% Point-by-point average of non-tweezer section 1

    if numNonTweezers1 > 0

        section1CurveMatrix = ...
            cell2mat(y(id,nonTweezer1Indices));

        section1ErrorMatrix = ...
            cell2mat(yerr(id,nonTweezer1Indices));

        [yNonTweezer1Mean{id}, ...
         yNonTweezer1Err{id}] = ...
            combineMeanMatrix( ...
                section1CurveMatrix, ...
                section1ErrorMatrix);

    else

        yNonTweezer1Mean{id} = nan(numX,1);
        yNonTweezer1Err{id} = nan(numX,1);
    end


    %% Point-by-point average of non-tweezer section 2

    if numNonTweezers2 > 0

        section2CurveMatrix = ...
            cell2mat(y(id,nonTweezer2Indices));

        section2ErrorMatrix = ...
            cell2mat(yerr(id,nonTweezer2Indices));

        [yNonTweezer2Mean{id}, ...
         yNonTweezer2Err{id}] = ...
            combineMeanMatrix( ...
                section2CurveMatrix, ...
                section2ErrorMatrix);

    else

        yNonTweezer2Mean{id} = nan(numX,1);
        yNonTweezer2Err{id} = nan(numX,1);
    end


    %% Scalar average across all real-tweezer ROI means

    [avgAllTweezers(id), ...
     avgAllTweezersErr(id)] = ...
        combineMeanValues( ...
            avgSignal(id,realTweezerIndices), ...
            avgSignalErr(id,realTweezerIndices));


    %% Scalar average across non-tweezer section 1 ROI means

    if numNonTweezers1 > 0

        [avgAllNonTweezers1(id), ...
         avgAllNonTweezers1Err(id)] = ...
            combineMeanValues( ...
                avgSignal(id,nonTweezer1Indices), ...
                avgSignalErr(id,nonTweezer1Indices));
    end


    %% Scalar average across non-tweezer section 2 ROI means

    if numNonTweezers2 > 0

        [avgAllNonTweezers2(id), ...
         avgAllNonTweezers2Err(id)] = ...
            combineMeanValues( ...
                avgSignal(id,nonTweezer2Indices), ...
                avgSignalErr(id,nonTweezer2Indices));
    end
end


%% Plot each individual ROI

scanColors = lines(numScanIDs);

for roiNum = 1:numROIs

    figure;
    hold on;

    for id = 1:numScanIDs

        if isempty(x{id})
            continue;
        end

        errorbar( ...
            x{id}, ...
            y{id,roiNum}, ...
            yerr{id,roiNum}, ...
            'o-', ...
            'Color',scanColors(id,:), ...
            'MarkerFaceColor',scanColors(id,:), ...
            'DisplayName',sprintf( ...
                '%s = %g', ...
                analyVar.avgScanParam, ...
                scanIDs(id)));
    end

    xlabel( ...
        analyVar.avgScanParam, ...
        'Interpreter','none');

    ylabel( ...
        depVarField, ...
        'Interpreter','none');

    title(sprintf('Tweezer ROI %d',roiNum));

    legend('Location','best');
    grid on;
    box on;
    hold off;
end


%% Plot average real-tweezer curve

plotGroupCurves( ...
    x, ...
    yTweezerMean, ...
    yTweezerErr, ...
    scanIDs, ...
    scanColors, ...
    analyVar.avgScanParam, ...
    depVarField, ...
    sprintf( ...
        'Mean %s across %d real tweezers', ...
        depVarField, ...
        numRealTweezers), ...
    'Average of Real-Tweezer ROI Signals');


%% Plot average non-tweezer section 1 curve

if numNonTweezers1 > 0

    plotGroupCurves( ...
        x, ...
        yNonTweezer1Mean, ...
        yNonTweezer1Err, ...
        scanIDs, ...
        scanColors, ...
        analyVar.avgScanParam, ...
        depVarField, ...
        sprintf( ...
            'Mean %s across %d section-1 ROIs', ...
            depVarField, ...
            numNonTweezers1), ...
        sprintf( ...
            'Average of Non-Tweezer Section 1: ROIs %d-%d', ...
            nonTweezer1Indices(1), ...
            nonTweezer1Indices(end)));
end


%% Plot average non-tweezer section 2 curve

if numNonTweezers2 > 0

    plotGroupCurves( ...
        x, ...
        yNonTweezer2Mean, ...
        yNonTweezer2Err, ...
        scanIDs, ...
        scanColors, ...
        analyVar.avgScanParam, ...
        depVarField, ...
        sprintf( ...
            'Mean %s across %d section-2 ROIs', ...
            depVarField, ...
            numNonTweezers2), ...
        sprintf( ...
            'Average of Non-Tweezer Section 2: ROIs %d-%d', ...
            nonTweezer2Indices(1), ...
            nonTweezer2Indices(end)));
end


%% Compare all three ROI groups for each scan ID

comparisonColors = lines(3);

for id = 1:numScanIDs

    if isempty(x{id})
        continue;
    end

    figure;
    hold on;

    errorbar( ...
        x{id}, ...
        yTweezerMean{id}, ...
        yTweezerErr{id}, ...
        'o-', ...
        'Color',comparisonColors(1,:), ...
        'MarkerFaceColor',comparisonColors(1,:), ...
        'DisplayName',sprintf( ...
            'Real tweezers: ROIs %d-%d', ...
            realTweezerIndices(1), ...
            realTweezerIndices(end)));

    if numNonTweezers1 > 0

        errorbar( ...
            x{id}, ...
            yNonTweezer1Mean{id}, ...
            yNonTweezer1Err{id}, ...
            's-', ...
            'Color',comparisonColors(2,:), ...
            'MarkerFaceColor',comparisonColors(2,:), ...
            'DisplayName',sprintf( ...
                'Non-tweezer section 1: ROIs %d-%d', ...
                nonTweezer1Indices(1), ...
                nonTweezer1Indices(end)));
    end

    if numNonTweezers2 > 0

        errorbar( ...
            x{id}, ...
            yNonTweezer2Mean{id}, ...
            yNonTweezer2Err{id}, ...
            'd-', ...
            'Color',comparisonColors(3,:), ...
            'MarkerFaceColor',comparisonColors(3,:), ...
            'DisplayName',sprintf( ...
                'Non-tweezer section 2: ROIs %d-%d', ...
                nonTweezer2Indices(1), ...
                nonTweezer2Indices(end)));
    end

    xlabel( ...
        analyVar.avgScanParam, ...
        'Interpreter','none');

    ylabel( ...
        depVarField, ...
        'Interpreter','none');

    title(sprintf( ...
        '%s = %g: ROI Group Comparison', ...
        analyVar.avgScanParam, ...
        scanIDs(id)), ...
        'Interpreter','none');

    legend('Location','best');
    grid on;
    box on;
    hold off;
end


%% Print scalar group averages

fprintf('\nScalar averages over each complete ROI curve:\n');

for id = 1:numScanIDs

    fprintf( ...
        ['Scan ID %g:\n' ...
         '  Real tweezers: %.6g +/- %.6g\n'], ...
        scanIDs(id), ...
        avgAllTweezers(id), ...
        avgAllTweezersErr(id));

    if numNonTweezers1 > 0

        fprintf( ...
            '  Non-tweezer section 1: %.6g +/- %.6g\n', ...
            avgAllNonTweezers1(id), ...
            avgAllNonTweezers1Err(id));
    end

    if numNonTweezers2 > 0

        fprintf( ...
            '  Non-tweezer section 2: %.6g +/- %.6g\n', ...
            avgAllNonTweezers2(id), ...
            avgAllNonTweezers2Err(id));
    end
end

end


function [combinedMean,combinedError] = ...
    combineMeanValues(values,errors)

% Combine scalar values using an equal-weight arithmetic mean.
% The uncertainty includes propagated input uncertainty and scatter
% between the input means.

values = values(:);
errors = errors(:);

valid = isfinite(values);

values = values(valid);
errors = errors(valid);

numValues = numel(values);

if numValues == 0

    combinedMean = NaN;
    combinedError = NaN;
    return;
end

combinedMean = mean(values);

% Missing individual errors cannot contribute to propagated uncertainty.
validErrors = isfinite(errors);

if any(validErrors)

    propagatedError = ...
        sqrt(sum(errors(validErrors).^2)) / numValues;

else

    propagatedError = 0;
end

if numValues > 1

    scatterError = ...
        std(values,0) / sqrt(numValues);

else

    scatterError = 0;
end

combinedError = sqrt( ...
    propagatedError.^2 + scatterError.^2);

end


function [meanCurve,errorCurve] = ...
    combineMeanMatrix(valueMatrix,errorMatrix)

% Combine ROI curves row-by-row using an equal-weight arithmetic mean.
% Each row corresponds to one independent-variable value.
% Each column corresponds to one ROI.

if ~isequal(size(valueMatrix),size(errorMatrix))

    error('The value and error matrices must have identical dimensions.');
end

numPoints = size(valueMatrix,1);

meanCurve = nan(numPoints,1);
errorCurve = nan(numPoints,1);

for pointNum = 1:numPoints

    [meanCurve(pointNum),errorCurve(pointNum)] = ...
        combineMeanValues( ...
            valueMatrix(pointNum,:), ...
            errorMatrix(pointNum,:));
end

end


function plotGroupCurves( ...
    x, ...
    meanCurves, ...
    errorCurves, ...
    scanIDs, ...
    scanColors, ...
    xLabelText, ...
    depVarField, ...
    yLabelText, ...
    titleText)

figure;
hold on;

for id = 1:numel(scanIDs)

    if isempty(x{id})
        continue;
    end

    errorbar( ...
        x{id}, ...
        meanCurves{id}, ...
        errorCurves{id}, ...
        'o-', ...
        'Color',scanColors(id,:), ...
        'MarkerFaceColor',scanColors(id,:), ...
        'DisplayName',sprintf( ...
            '%s = %g', ...
            xLabelText, ...
            scanIDs(id)));
end

xlabel( ...
    xLabelText, ...
    'Interpreter','none');

ylabel( ...
    yLabelText, ...
    'Interpreter','none');

title( ...
    titleText, ...
    'Interpreter','none');

legend('Location','best');
grid on;
box on;
hold off;

end