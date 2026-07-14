function avgDataset = add_fit_avg_batch(analyVar, avgDataset)
% Adds averaged OD images and averaged fit parameters into avgDataset.
%
% Adds:
%   avgDataset.All_OD_Image
%   avgDataset.All_PCell
%   avgDataset.All_fitParams

%% Preallocate
[avgDataset.All_OD_Image, ...
 avgDataset.All_PCell, ...
 avgDataset.All_fitParams] = deal(cell(1, avgDataset.CounterAtom));

%% Rebuild paramFitFiles if missing
if ~isfield(avgDataset,'paramFitFiles')
    avgDataset.paramFitFiles = cell(avgDataset.CounterAtom,1);

    for j = 1:avgDataset.CounterAtom
        safeVal = regexprep(num2str(avgDataset.imagevcoAtom(ceil(j/avgDataset.numTweezers)),'%.12g'), ...
            '[^a-zA-Z0-9_\-\.]', '_');

        avgDataset.paramFitFiles{j} = fullfile(avgDataset.avgDir, ...
            ['AvgFit_' analyVar.avgScanParamField '_' safeVal '_' ...
            analyVar.sampleType analyVar.fitModel ...
            analyVar.paramFitFilename analyVar.paramFitFileExt]);
    end
end

%% SubPlot Info (sizing)
[avgDataset.SubPlotRows, avgDataset.SubPlotCols] = optiSubPlotNum(avgDataset.CounterAtom);

%% Loop through averaged images
for j = 1:avgDataset.CounterAtom

    %% Read averaged OD image
    if exist(avgDataset.avgODFiles{j}, 'file')
        avgDataset.All_OD_Image{j} = dlmread(avgDataset.avgODFiles{j});
    else
        error('imagefit:NoAvgODSaved', ...
            'Cannot load averaged OD image:\n%s', avgDataset.avgODFiles{j});
    end

    %% Read averaged fit parameters
    fitFile = avgDataset.paramFitFiles{j};

    if exist(fitFile, 'file')

        rawFit = dlmread(fitFile);

        % Same logic as original:
        % rows = fitted lattice axes/cloud regions
        % cols = fit params + weighting flag
        avgDataset.All_PCell{j} = mat2cell( ...
            rawFit, ...
            ones(1, sum(analyVar.LatticeAxesFit)), ...
            length(analyVar.InitCondList) + 1);

        % Convert numerical fit params into struct fields
        avgDataset.All_fitParams{j} = cellfun( ...
            @(P) cell2struct( ...
                num2cell(P(1:length(analyVar.InitCondList))), ...
                analyVar.InitCondList, ...
                2), ...
            avgDataset.All_PCell{j}, ...
            'UniformOutput', 0);

    else
        error('imagefit:NoAvgFitSaved', ...
            'Cannot load averaged fit file:\n%s', fitFile);
    end
end

%% Save updated avgDataset
save(fullfile(avgDataset.avgDir, 'avgDataset.mat'), 'avgDataset');

end