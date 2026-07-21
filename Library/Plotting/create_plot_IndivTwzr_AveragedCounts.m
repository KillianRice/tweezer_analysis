function create_plot_IndivTwzr_AveragedCounts(analyVar, indivDataset, avgDataset, plotRawCounts)
    %%% average_plot.m - Joe Whalen 2017.12.15
    %%% Make a plot of the average of any quantity in indivDataset grouped
    %%% by the flags in the master batch file.
    
    indVarField = 'imagevcoAtom'; % The Field of an IndivDataset that is to be plotted on the X axis
    if plotRawCounts
        depVarField = 'OD_TotalCountsImg1Raw';
    else
        depVarField = 'OD_TotalCounts';  % The field of an indivdataset that is to be plotted on the y axis
    end
    
    [xdata_clean, ydata_clean] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    scanIDs = analyVar.uniqScanList;
    
    %% Determine number of tweezers
    load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
    numTweezers = size(tweezerROI.centersXY,1);
    
    %% Preallocate
    numScanIDs = numel(scanIDs);
    
    x            = cell(numScanIDs,1);
    y            = cell(numScanIDs,numTweezers);
    yerr         = cell(numScanIDs,numTweezers);
    
    AvgSig       = nan(numScanIDs,numTweezers);
    AvgSig_err   = nan(numScanIDs,numTweezers);
    
    yTweezerMean = cell(numScanIDs,1);
    yTweezerErr  = cell(numScanIDs,1);

    ynonTweezerMean = cell(numScanIDs,1);
    ynonTweezerErr  = cell(numScanIDs,1);
    
    AvgAllTweezers     = nan(numScanIDs,1);
    AvgAllTweezers_err = nan(numScanIDs,1);
    
    AvgAllNonTweezers     = nan(numScanIDs,1);
    AvgAllNonTweezers_err = nan(numScanIDs,1);

    fprintf('Plotting Individual Tweezer Counts.\n\n')
    
    %% Average repeated scans separately for each tweezer
    for id = 1:numScanIDs
    
        %% Collect all x values belonging to this scan ID
        xAll = [];
    
        for basename = 1:analyVar.numBasenamesAtom
    
            if scanIDs(id) ~= analyVar.meanListVar(basename)
                continue;
            end
    
            xBase = xdata_clean{basename};
            xBase = xBase(:);
    
            xAll = [xAll; xBase];
        end
    
        % Unique sorted independent-variable values
        x{id} = unique(xAll, 'sorted');
    
        numX = numel(x{id});
    
        for tweezerNum = 1:numTweezers
            y{id,tweezerNum}    = nan(numX,1);
            yerr{id,tweezerNum} = nan(numX,1);
        end
    
        %% Loop through each independent-variable value
        for i = 1:numX
    
            currentX = x{id}(i);
    
            % One collection vector per tweezer
            valuesAtX = cell(1,numTweezers);
    
            %% Search every matching input file
            for basename = 1:analyVar.numBasenamesAtom
    
                if scanIDs(id) ~= analyVar.meanListVar(basename)
                    continue;
                end
    
                xBase = xdata_clean{basename}(:);
                yBase = ydata_clean{basename};
    
                if size(yBase,1) ~= numel(xBase)
                    error(['Number of rows in ydata_clean{%d} does not match ' ...
                           'the number of x values.'], basename);
                end
    
                % Tolerance-based matching is safer than exact floating equality
                xTol = max(1e-12, 1e-9*max(1,abs(currentX)));
                matchingRows = find(abs(xBase-currentX) <= xTol);
    
                for rowNum = matchingRows(:)'
    
                    for tweezerNum = 1:min(numTweezers,size(yBase,2))
    
                        thisValue = yBase(rowNum,tweezerNum);
    
                        if ~isnan(thisValue)
                            valuesAtX{tweezerNum}(end+1,1) = thisValue;
                        end
                    end
                end
            end
    
            %% Average repeated-file values for each tweezer
            for tweezerNum = 1:numTweezers
    
                vals = valuesAtX{tweezerNum};
    
                if isempty(vals)
                    continue;
                end
    
                y{id,tweezerNum}(i) = mean(vals,'omitnan');
    
                nValid = sum(~isnan(vals));
    
                if nValid > 1
                    yerr{id,tweezerNum}(i) = ...
                        std(vals,'omitnan') / sqrt(nValid);
                elseif nValid == 1
                    yerr{id,tweezerNum}(i) = 0;
                end
            end
        end
    
        %% Scalar average across the independent variable for each tweezer
        for tweezerNum = 1:numTweezers
    
            thisCurve = y{id,tweezerNum};
    
            AvgSig(id,tweezerNum) = mean(thisCurve,'omitnan');
    
            nValid = sum(~isnan(thisCurve));
    
            if nValid > 1
                AvgSig_err(id,tweezerNum) = ...
                    std(thisCurve,'omitnan') / sqrt(nValid);
            elseif nValid == 1
                AvgSig_err(id,tweezerNum) = 0;
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define ROI groups
    
    numNonTweezers1 = analyVar.numFakeTweezers1;
    numNonTweezers2 = analyVar.numFakeTweezers2;
    
    numRealTweezers = ...
        numTweezers - numNonTweezers1 - numNonTweezers2;
    
    %% Validate the requested grouping
    
    groupCounts = [
        numRealTweezers
        numNonTweezers1
        numNonTweezers2
    ];
    
    if any(groupCounts < 0) || any(groupCounts ~= round(groupCounts))
        error(['The numbers of real tweezers and non-tweezer ROIs ' ...
               'must all be nonnegative integers.']);
    end
    
    if numRealTweezers < 1
        error(['The selected non-tweezer section sizes leave no real ' ...
               'tweezer ROIs.']);
    end
    
    if sum(groupCounts) ~= numTweezers
        error(['ROI grouping error: %d real + %d section-1 + ' ...
               '%d section-2 does not equal %d total ROIs.'], ...
            numRealTweezers, ...
            numNonTweezers1, ...
            numNonTweezers2, ...
            numTweezers);
    end
    
    %% ROI indices
    %
    % Example for 12 total, section 1 = 2, section 2 = 4:
    %   realTweezerIndices = 1:6
    %   nonTweezer1Indices = 7:8
    %   nonTweezer2Indices = 9:12
    
    realTweezerIndices = ...
        1:numRealTweezers;
    
    nonTweezer1Indices = ...
        (numRealTweezers + 1): ...
        (numRealTweezers + numNonTweezers1);
    
    nonTweezer2Indices = ...
        (numRealTweezers + numNonTweezers1 + 1): ...
        numTweezers;
    
    
    %% Average the real-tweezer curves
    
    tweezerCurveMatrix = ...
        nan(numX,numRealTweezers);
    
    for localIndex = 1:numRealTweezers
    
        roiNum = realTweezerIndices(localIndex);
    
        tweezerCurveMatrix(:,localIndex) = ...
            y{id,roiNum};
    end
    
    yTweezerMean{id} = ...
        mean(tweezerCurveMatrix,2,'omitnan');
    
    nTweezersAtPoint = ...
        sum(~isnan(tweezerCurveMatrix),2);
    
    yTweezerErr{id} = ...
        std(tweezerCurveMatrix,0,2,'omitnan') ./ ...
        sqrt(nTweezersAtPoint);
    
    yTweezerErr{id}(nTweezersAtPoint <= 1) = 0;
    
    
    %% Average non-tweezer section 1
    
    if numNonTweezers1 > 0
    
        nonTweezer1CurveMatrix = ...
            nan(numX,numNonTweezers1);
    
        for localIndex = 1:numNonTweezers1
    
            roiNum = nonTweezer1Indices(localIndex);
    
            nonTweezer1CurveMatrix(:,localIndex) = ...
                y{id,roiNum};
        end
    
        yNonTweezer1Mean{id} = ...
            mean(nonTweezer1CurveMatrix,2,'omitnan');
    
        nNonTweezers1AtPoint = ...
            sum(~isnan(nonTweezer1CurveMatrix),2);
    
        yNonTweezer1Err{id} = ...
            std(nonTweezer1CurveMatrix,0,2,'omitnan') ./ ...
            sqrt(nNonTweezers1AtPoint);
    
        yNonTweezer1Err{id}(nNonTweezers1AtPoint <= 1) = 0;
    
    else
    
        yNonTweezer1Mean{id} = nan(numX,1);
        yNonTweezer1Err{id} = nan(numX,1);
    end
    
    
    %% Average non-tweezer section 2
    
    if numNonTweezers2 > 0
    
        nonTweezer2CurveMatrix = ...
            nan(numX,numNonTweezers2);
    
        for localIndex = 1:numNonTweezers2
    
            roiNum = nonTweezer2Indices(localIndex);
    
            nonTweezer2CurveMatrix(:,localIndex) = ...
                y{id,roiNum};
        end
    
        yNonTweezer2Mean{id} = ...
            mean(nonTweezer2CurveMatrix,2,'omitnan');
    
        nNonTweezers2AtPoint = ...
            sum(~isnan(nonTweezer2CurveMatrix),2);
    
        yNonTweezer2Err{id} = ...
            std(nonTweezer2CurveMatrix,0,2,'omitnan') ./ ...
            sqrt(nNonTweezers2AtPoint);
    
        yNonTweezer2Err{id}(nNonTweezers2AtPoint <= 1) = 0;
    
    else
    
        yNonTweezer2Mean{id} = nan(numX,1);
        yNonTweezer2Err{id} = nan(numX,1);
    end
    
    
    %% One scalar per scan ID: real tweezers
    
    realTweezerAvgValues = ...
        AvgSig(id,realTweezerIndices);
    
    AvgAllTweezers(id) = ...
        mean(realTweezerAvgValues,'omitnan');
    
    nValidTweezers(id) = ...
        sum(~isnan(realTweezerAvgValues));
    
    if nValidTweezers(id) > 1
    
        AvgAllTweezers_err(id) = ...
            std(realTweezerAvgValues,0,'omitnan') / ...
            sqrt(nValidTweezers(id));
    
    elseif nValidTweezers(id) == 1
    
        AvgAllTweezers_err(id) = 0;
    
    else
    
        AvgAllTweezers_err(id) = NaN;
    end
    
    
    %% One scalar per scan ID: non-tweezer section 1
    
    if numNonTweezers1 > 0
    
        nonTweezer1AvgValues = ...
            AvgSig(id,nonTweezer1Indices);
    
        AvgAllNonTweezers1(id) = ...
            mean(nonTweezer1AvgValues,'omitnan');
    
        nValidNonTweezers1(id) = ...
            sum(~isnan(nonTweezer1AvgValues));
    
        if nValidNonTweezers1(id) > 1
    
            AvgAllNonTweezers1_err(id) = ...
                std(nonTweezer1AvgValues,0,'omitnan') / ...
                sqrt(nValidNonTweezers1(id));
    
        elseif nValidNonTweezers1(id) == 1
    
            AvgAllNonTweezers1_err(id) = 0;
    
        else
    
            AvgAllNonTweezers1_err(id) = NaN;
        end
    
    else
    
        AvgAllNonTweezers1(id) = NaN;
        AvgAllNonTweezers1_err(id) = NaN;
        nValidNonTweezers1(id) = 0;
    end
    
    
    %% One scalar per scan ID: non-tweezer section 2
    
    if numNonTweezers2 > 0
    
        nonTweezer2AvgValues = ...
            AvgSig(id,nonTweezer2Indices);
    
        AvgAllNonTweezers2(id) = ...
            mean(nonTweezer2AvgValues,'omitnan');
    
        nValidNonTweezers2(id) = ...
            sum(~isnan(nonTweezer2AvgValues));
    
        if nValidNonTweezers2(id) > 1
    
            AvgAllNonTweezers2_err(id) = ...
                std(nonTweezer2AvgValues,0,'omitnan') / ...
                sqrt(nValidNonTweezers2(id));
    
        elseif nValidNonTweezers2(id) == 1
    
            AvgAllNonTweezers2_err(id) = 0;
    
        else
    
            AvgAllNonTweezers2_err(id) = NaN;
        end
    
    else
    
        AvgAllNonTweezers2(id) = NaN;
        AvgAllNonTweezers2_err(id) = NaN;
        nValidNonTweezers2(id) = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot each individual Averaged Tweezer
    colors = lines(numScanIDs);

    for tweezerNum = 1:numTweezers
    
        figure;
        hold on;
    
        legendList = cell(numScanIDs,1);
    
        for id = 1:numScanIDs
    
            errorbar( ...
                x{id}, ...
                y{id,tweezerNum}, ...
                yerr{id,tweezerNum}, ...
                'o-', ...
                'Color', colors(id,:), ...
                'MarkerFaceColor', colors(id,:));
    
            legendList{id} = sprintf('%s = %g', ...
                analyVar.avgScanParam, scanIDs(id));
        end
    
        xlabel(analyVar.avgScanParam, 'Interpreter','none');
        ylabel(depVarField, 'Interpreter','none');
    
        title(sprintf('Tweezer ROI %d', tweezerNum));
    
        legend(legendList,'Location','best');
        grid on;
        hold off;
    end
    
    %% Plot the average of the real-tweezer spots
    
    figure;
    hold on;
    
    legendList = cell(numScanIDs,1);
    colors = lines(numScanIDs);
    
    for id = 1:numScanIDs
    
        errorbar( ...
            x{id}, ...
            yTweezerMean{id}, ...
            yTweezerErr{id}, ...
            'o-', ...
            'Color',colors(id,:), ...
            'MarkerFaceColor',colors(id,:));
    
        legendList{id} = sprintf('%s = %g', ...
            analyVar.avgScanParam,scanIDs(id));
    end
    
    xlabel(analyVar.avgScanParam,'Interpreter','none');
    
    ylabel(sprintf( ...
        'Mean %s across %d real tweezers', ...
        depVarField,numRealTweezers), ...
        'Interpreter','none');
    
    title('Average of Real-Tweezer ROI Signals');
    
    legend(legendList,'Location','best');
    grid on;
    hold off;
    
    
    %% Plot the average of non-tweezer section 1
    
    if numNonTweezers1 > 0
    
        figure;
        hold on;
    
        legendList = cell(numScanIDs,1);
    
        for id = 1:numScanIDs
    
            errorbar( ...
                x{id}, ...
                yNonTweezer1Mean{id}, ...
                yNonTweezer1Err{id}, ...
                'o-', ...
                'Color',colors(id,:), ...
                'MarkerFaceColor',colors(id,:));
    
            legendList{id} = sprintf('%s = %g', ...
                analyVar.avgScanParam,scanIDs(id));
        end
    
        xlabel(analyVar.avgScanParam,'Interpreter','none');
    
        ylabel(sprintf( ...
            'Mean %s across %d section-1 ROIs', ...
            depVarField,numNonTweezers1), ...
            'Interpreter','none');
    
        title(sprintf( ...
            'Average of Non-Tweezer Section 1: ROIs %d-%d', ...
            nonTweezer1Indices(1),nonTweezer1Indices(end)));
    
        legend(legendList,'Location','best');
        grid on;
        hold off;
    end
    
    
    %% Plot the average of non-tweezer section 2
    
    if numNonTweezers2 > 0
    
        figure;
        hold on;
    
        legendList = cell(numScanIDs,1);
    
        for id = 1:numScanIDs
    
            errorbar( ...
                x{id}, ...
                yNonTweezer2Mean{id}, ...
                yNonTweezer2Err{id}, ...
                'o-', ...
                'Color',colors(id,:), ...
                'MarkerFaceColor',colors(id,:));
    
            legendList{id} = sprintf('%s = %g', ...
                analyVar.avgScanParam,scanIDs(id));
        end
    
        xlabel(analyVar.avgScanParam,'Interpreter','none');
    
        ylabel(sprintf( ...
            'Mean %s across %d section-2 ROIs', ...
            depVarField,numNonTweezers2), ...
            'Interpreter','none');
    
        title(sprintf( ...
            'Average of Non-Tweezer Section 2: ROIs %d-%d', ...
            nonTweezer2Indices(1),nonTweezer2Indices(end)));
    
        legend(legendList,'Location','best');
        grid on;
        hold off;
    end
    
    
    %% Compare all three ROI-group averages for each scan ID
    
    comparisonColors = lines(3);
    
    for id = 1:numScanIDs
    
        figure;
        hold on;
    
        comparisonLegend = {};
    
        %% Real-tweezer group
    
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
    
        comparisonLegend{end+1} = sprintf( ...
            'Real tweezers: ROIs %d-%d', ...
            realTweezerIndices(1), ...
            realTweezerIndices(end));
    
    
        %% Non-tweezer section 1
    
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
    
            comparisonLegend{end+1} = sprintf( ...
                'Non-tweezer section 1: ROIs %d-%d', ...
                nonTweezer1Indices(1), ...
                nonTweezer1Indices(end));
        end
    
    
        %% Non-tweezer section 2
    
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
    
            comparisonLegend{end+1} = sprintf( ...
                'Non-tweezer section 2: ROIs %d-%d', ...
                nonTweezer2Indices(1), ...
                nonTweezer2Indices(end));
        end
    
    
        xlabel(analyVar.avgScanParam,'Interpreter','none');
        ylabel(depVarField,'Interpreter','none');
    
        title(sprintf( ...
            '%s = %g: ROI Group Comparison', ...
            analyVar.avgScanParam,scanIDs(id)), ...
            'Interpreter','none');
    
        legend(comparisonLegend,'Location','best');
    
        grid on;
        box on;
        hold off;
    end
    
end