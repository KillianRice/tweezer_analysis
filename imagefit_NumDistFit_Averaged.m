function imagefit_NumDistFit_Averaged(varargin)

%% Load variables and file data
if nargin == 0
    analyVar = AnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
    load(fullfile(analyVar.analyOutDir, analyVar.avgOutSubDir, 'avgDataset.mat'),'avgDataset');
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
    avgDataset = varargin{3};
end

avgDataset.paramFitFiles = cell(avgDataset.CounterAtom,1);

fprintf('\nFitting averaged OD images...\n'); 

for j = 1:avgDataset.CounterAtom

    fprintf('Fitting averaged image %d of %d\n', j, avgDataset.CounterAtom);

    %%Grab indivDataset info within avgDataset
    src = avgDataset.sourceInfo{j};

    basenameNum = src.basenameNums(1);
    k           = src.imageNums(1);
    tweezerNum  = src.tweezerNums(1);

    OD_Image_Single = dlmread(avgDataset.avgODFiles{j});

    if analyVar.fitSmoothOD
        OD_Image_Single = analyVar.smoothFilt(OD_Image_Single, analyVar.smoothFiltMat);
    end

    % Use first original dataset as geometry reference
    % basenameNum = 1;
    %roiWin_Index = indivDataset{1}.roiWin_Index;%get_avg_roiWin_Index(analyVar, OD_Image_Single);
    if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
    
        roiWin_Index = {true(size(OD_Image_Single))};
        fitWinSize = size(OD_Image_Single,1);

        OD_Fit_ImageCell = cellfun( ...
            @(x) reshape(OD_Image_Single(x), [1 1]*fitWinSize), ...
            roiWin_Index, ...
            'UniformOutput', 0);
    
    else
    
        roiWin_Index = avgDataset.roiWin_Index;
        fitWinSize = analyVar.funcFitWin(basenameNum);

        OD_Fit_ImageCell = cellfun( ...
            @(x) reshape(OD_Image_Single(x), [1 1]*analyVar.funcFitWin(basenameNum)), ...
            cellfun(analyVar.fitWinLogicInd, roiWin_Index, 'UniformOutput', 0), ...
            'UniformOutput', 0);
    
    end
    
    % For averaged images, j is replaced by k
    InitGuess = get_fit_params(analyVar, OD_Image_Single, OD_Fit_ImageCell, indivDataset, basenameNum, k);

    [lowBndVec, upBndVec] = get_fit_bounds(analyVar, basenameNum);

    nFits = numel(InitGuess);
    
    lowBndCell = cell(size(InitGuess));
    upBndCell  = cell(size(InitGuess));
    
    for q = 1:nFits
        lowBndCell{q} = lowBndVec(:);
        upBndCell{q}  = upBndVec(:);
    end

    errCell = cellfun(@(x) get_OD_weight(analyVar.weightPeak,x), ...
        OD_Fit_ImageCell, 'UniformOutput', 0);

    %[Xgrid,Ygrid] = meshgrid(1:length(OD_Fit_ImageCell{1}));
    [Xgrid,Ygrid] = meshgrid(1:fitWinSize);

    optimOpt = optimset( ...
        'Display','off', ...
        'FinDiffType','central', ...
        'TolFun',1e-9, ...
        'TolX',1e-9);

    PCell = cellfun( ...
        @(x,y,z,lb,ub) lsqcurvefit( ...
            str2func(analyVar.fitModel), ...
            z, ...
            [Xgrid(:), Ygrid(:), y(:)], ...
            x(:), ...
            lb, ...
            ub, ...
            optimOpt), ...
        cellfun(@(m,n) m(:)./n(:), OD_Fit_ImageCell, errCell, 'UniformOutput', 0), ...
        errCell, ...
        InitGuess, ...
        lowBndCell, ...
        upBndCell, ...
        'UniformOutput', 0);

    PCell = cellfun(@(x) [x' analyVar.weightPeak], PCell, 'UniformOutput', 0);
    
    %%% Use ceil(j/avgDataset.numTweezers) since each tweezer within the
    %%% same image should correlate to the same independent variable
    safeVal = regexprep(num2str(avgDataset.imagevcoAtom(k),'%.12g'),'[^a-zA-Z0-9_\-\.]','_');

    paramFile = fullfile(avgDataset.avgDir, ...
        ['AvgFit_' analyVar.avgScanParamField '_' safeVal '_' ...
        analyVar.sampleType analyVar.fitModel analyVar.paramFitFilename ...
        analyVar.paramFitFileExt]);

    dlmwrite(paramFile, PCell, '\t');

    avgDataset.paramFitFiles{j} = paramFile;

end

save(fullfile(avgDataset.avgDir,'avgDataset.mat'),'avgDataset');

fclose('all');

fprintf('Averaged cloud fitting completed.\n\n');

end