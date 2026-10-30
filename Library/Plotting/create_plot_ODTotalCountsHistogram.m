function create_plot_ODTotalCountsHistogram(analyVar, avgDataset, plotIndividual)

plotIndividual = 1;

%% Group IDs and Scam IDs
groupVals = avgDataset.imagevcoAtom;
entryScanID = avgDataset.scanID;

%% Identify unique group ID and scan ID combinations

groupScanVals = unique([groupVals, entryScanID], 'rows');
uniqueGroups = groupScanVals(:,1);
uniqueScanIDs = groupScanVals(:,2);

%% Bin edges

allCounts = [];

for j = 1:avgDataset.CounterAtom
    if analyVar.plotRawCounts
        allCounts = [allCounts; avgDataset.sourceInfo{j}.rawODCounts(:)];
    else
        allCounts = [allCounts; avgDataset.sourceInfo{j}.totODCounts(:)];
    end
end

edges = linspace(min(allCounts), max(allCounts), 100);

% Combined plot
figure;
hold on;

colors = lines(size(groupScanVals,1));
legendList = cell(size(groupScanVals,1),1);

for g = 1:size(groupScanVals,1)

    groupID = uniqueGroups(g);
    scanID = uniqueScanIDs(g);
    counts = [];

    for j = 1:avgDataset.CounterAtom
        if groupVals(j) == groupID && entryScanID(j) == scanID
            if analyVar.plotRawCounts
                counts = [counts; avgDataset.IndivTwzrRawCounts{j}(:)];
            else
                counts = [counts; avgDataset.IndivTwzrTotODCounts{j}(:)];
            end
        end
    end

    histogram(counts, ...
        'BinEdges', edges, ...
        'Normalization', 'probability', ...
        'FaceColor', colors(g,:), ...
        'FaceAlpha', 0.45);

    legendList{g} = sprintf('%s = %g, Scan ID = %g', ...
        analyVar.avgScanParam, groupID, scanID);
end

xlabel('Integrated OD counts');
ylabel('Probability');
title('OD Total Counts Histogram by Group and Scan ID');
grid on;
legend(legendList, 'Location', 'best');
hold off;


%% Optional individual plots
if plotIndividual

    for g = 1:numel(uniqueGroups)

        groupID = uniqueGroups(g);

        % Find all scan IDs associated with this group ID
        scanIDsForGroup = uniqueScanIDs(uniqueGroups == groupID);

        figure;
        hold on;

        groupColors = lines(numel(scanIDsForGroup));
        legendList = cell(numel(scanIDsForGroup),1);

        for s = 1:numel(scanIDsForGroup)

            scanID = scanIDsForGroup(s);
            counts = [];

            for j = 1:avgDataset.CounterAtom
                if groupVals(j) == groupID && entryScanID(j) == scanID
                    if analyVar.plotRawCounts
                        counts = [counts; avgDataset.IndivTwzrRawCounts{j}(:)];
                    else
                        counts = [counts; avgDataset.IndivTwzrTotODCounts{j}(:)];
                    end
                end
            end

            histogram(counts, ...
                'BinEdges', edges, ...
                'Normalization', 'probability', ...
                'FaceColor', groupColors(s,:), ...
                'FaceAlpha', 0.45);

            legendList{s} = sprintf('Scan ID = %g', scanID);

        end

        xlabel('Integrated OD counts');
        ylabel('Probability');
        title(sprintf('OD Total Counts Histogram - %s = %g', ...
            analyVar.avgScanParam, groupID));
        grid on;
        legend(legendList, 'Location', 'best');
        hold off;

    end
end