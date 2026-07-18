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
    
    AvgAllTweezers     = nan(numScanIDs,1);
    AvgAllTweezers_err = nan(numScanIDs,1);

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
    
        %% Average the already-averaged tweezer curves
        tweezerCurveMatrix = nan(numX,numTweezers);
    
        for tweezerNum = 1:numTweezers
            tweezerCurveMatrix(:,tweezerNum) = y{id,tweezerNum};
        end
    
        yTweezerMean{id} = mean(tweezerCurveMatrix,2,'omitnan');
    
        nTweezersAtPoint = sum(~isnan(tweezerCurveMatrix),2);
    
        yTweezerErr{id} = ...
            std(tweezerCurveMatrix,0,2,'omitnan') ./ sqrt(nTweezersAtPoint);
    
        yTweezerErr{id}(nTweezersAtPoint <= 1) = 0;
    
        %% One scalar per scan ID after combining all tweezers
        AvgAllTweezers(id) = mean(AvgSig(id,:),'omitnan');
    
        nValidTweezers = sum(~isnan(AvgSig(id,:)));
    
        if nValidTweezers > 1
            AvgAllTweezers_err(id) = ...
                std(AvgSig(id,:),0,2,'omitnan') / sqrt(nValidTweezers);
        elseif nValidTweezers == 1
            AvgAllTweezers_err(id) = 0;
        end
    end
    

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
    
        xlabel(indVarField, 'Interpreter','none');
        ylabel(depVarField, 'Interpreter','none');
    
        title(sprintf('Tweezer ROI %d', tweezerNum));
    
        legend(legendList,'Location','best');
        grid on;
        hold off;
    end
    
    %% Plot the average of the Tweezer Spots
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
            'Color', colors(id,:), ...
            'MarkerFaceColor', colors(id,:));
    
        legendList{id} = sprintf('%s = %g', ...
            analyVar.avgScanParam, scanIDs(id));
    end
    
    xlabel(indVarField, 'Interpreter','none');
    ylabel(sprintf('Mean %s across tweezers', depVarField), ...
        'Interpreter','none');
    
    title('Average of Individually Averaged Tweezer Signals');
    
    legend(legendList,'Location','best');
    grid on;
    hold off;

    
end