function create_plot_ODTotalCountsHistogram(analyVar, avgDataset, plotIndividual)

if nargin < 3
    plotIndividual = 0;
end

%% Collect group IDs
groupVals = zeros(avgDataset.CounterAtom,1);

for j = 1:avgDataset.CounterAtom
    groupVals(j) = avgDataset.sourceInfo{j}.imagevcoAtom;
end

uniqueGroups = unique(groupVals);

%% Bin edges
if isfield(analyVar,'ODCountHistEdges')
    edges = analyVar.ODCountHistEdges;
else
    allCounts = [];

    for j = 1:avgDataset.CounterAtom
        allCounts = [allCounts; avgDataset.sourceInfo{j}.totODCounts(:)];
    end

    edges = linspace(min(allCounts), max(allCounts), 40);
end

%% Combined plot
figure;
hold on;

colors = lines(numel(uniqueGroups));
legendList = cell(numel(uniqueGroups),1);

for g = 1:numel(uniqueGroups)

    groupID = uniqueGroups(g);
    counts = [];

    for j = 1:avgDataset.CounterAtom
        if groupVals(j) == groupID
            counts = [counts; avgDataset.sourceInfo{j}.totODCounts(:)];
        end
    end

    histogram(counts, ...
        'BinEdges', edges, ...
        'Normalization', 'probability', ...
        'FaceColor', colors(g,:), ...
        'FaceAlpha', 0.45);

    legendList{g} = sprintf('%s = %g', analyVar.avgScanParam, groupID);
end

xlabel('Integrated OD counts');
ylabel('Probability');
title('OD Total Counts Histogram by avgDataset Group');
grid on;
legend(legendList, 'Location', 'best');
hold off;

%% Optional individual plots
if plotIndividual

    for g = 1:numel(uniqueGroups)

        groupID = uniqueGroups(g);
        counts = [];

        for j = 1:avgDataset.CounterAtom
            if groupVals(j) == groupID
                counts = [counts; avgDataset.sourceInfo{j}.totODCounts(:)];
            end
        end

        figure;
        histogram(counts, ...
            'BinEdges', edges, ...
            'Normalization', 'probability', ...
            'FaceColor', colors(g,:), ...
            'FaceAlpha', 0.6);

        xlabel('Integrated OD counts');
        ylabel('Probability');
        title(sprintf('OD Total Counts Histogram - %s = %g', ...
            analyVar.avgScanParam, groupID));
        grid on;

    end
end

end