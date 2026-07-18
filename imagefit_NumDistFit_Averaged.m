% function imagefit_NumDistFit_Averaged(varargin)
% 
% %% Load variables and file data
% if nargin == 0
%     analyVar = AnalysisVariables;
%     indivDataset = get_indiv_batch_data(analyVar);
%     load(fullfile(analyVar.analyOutDir, analyVar.avgOutSubDir, 'avgDataset.mat'),'avgDataset');
% else
%     analyVar     = varargin{1}; % if arguments are passed analyVar must be first
%     indivDataset = varargin{2}; % indivDataset must be second
%     avgDataset = varargin{3};
% end
% 
% avgDataset.paramFitFiles = cell(avgDataset.CounterAtom,1);
% 
% fprintf('\nFitting averaged OD images...\n'); 
% 
% for j = 1:avgDataset.CounterAtom
% 
%     fprintf('Fitting averaged image %d of %d\n', j, avgDataset.CounterAtom);
% 
%     %%Grab indivDataset info within avgDataset
%     src = avgDataset.sourceInfo{j};
% 
%     basenameNum = src.basenameNums(1);
%     k           = src.imageNums(1);
%     tweezerNum  = src.tweezerNums(1);
% 
%     OD_Image_Single = dlmread(avgDataset.avgODFiles{j});
% 
%     if analyVar.fitSmoothOD
%         OD_Image_Single = analyVar.smoothFilt(OD_Image_Single, analyVar.smoothFiltMat);
%     end
% 
%     % Use first original dataset as geometry reference
%     % basenameNum = 1;
%     %roiWin_Index = indivDataset{1}.roiWin_Index;%get_avg_roiWin_Index(analyVar, OD_Image_Single);
%     if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
%     
%         roiWin_Index = {true(size(OD_Image_Single))};
%         fitWinSize = size(OD_Image_Single,1);
% 
%         OD_Fit_ImageCell = cellfun( ...
%             @(x) reshape(OD_Image_Single(x), [1 1]*fitWinSize), ...
%             roiWin_Index, ...
%             'UniformOutput', 0);
%     
%     else
%     
%         roiWin_Index = avgDataset.roiWin_Index;
%         fitWinSize = analyVar.funcFitWin(basenameNum);
% 
%         OD_Fit_ImageCell = cellfun( ...
%             @(x) reshape(OD_Image_Single(x), [1 1]*analyVar.funcFitWin(basenameNum)), ...
%             cellfun(analyVar.fitWinLogicInd, roiWin_Index, 'UniformOutput', 0), ...
%             'UniformOutput', 0);
%     
%     end
%     
%     % For averaged images, j is replaced by k
%     InitGuess = get_fit_params(analyVar, OD_Image_Single, OD_Fit_ImageCell, indivDataset, basenameNum, k);
% 
%     [lowBndVec, upBndVec] = get_fit_bounds(analyVar, basenameNum);
% 
%     nFits = numel(InitGuess);
%     
%     lowBndCell = cell(size(InitGuess));
%     upBndCell  = cell(size(InitGuess));
%     
%     for q = 1:nFits
%         lowBndCell{q} = lowBndVec(:);
%         upBndCell{q}  = upBndVec(:);
%     end
% 
%     errCell = cellfun(@(x) get_OD_weight(analyVar.weightPeak,x), ...
%         OD_Fit_ImageCell, 'UniformOutput', 0);
% 
%     %[Xgrid,Ygrid] = meshgrid(1:length(OD_Fit_ImageCell{1}));
%     [Xgrid,Ygrid] = meshgrid(1:fitWinSize);
% 
%     optimOpt = optimset( ...
%         'Display','off', ...
%         'FinDiffType','central', ...
%         'TolFun',1e-9, ...
%         'TolX',1e-9);
% 
%     PCell = cellfun( ...
%         @(x,y,z,lb,ub) lsqcurvefit( ...
%             str2func(analyVar.fitModel), ...
%             z, ...
%             [Xgrid(:), Ygrid(:), y(:)], ...
%             x(:), ...
%             lb, ...
%             ub, ...
%             optimOpt), ...
%         cellfun(@(m,n) m(:)./n(:), OD_Fit_ImageCell, errCell, 'UniformOutput', 0), ...
%         errCell, ...
%         InitGuess, ...
%         lowBndCell, ...
%         upBndCell, ...
%         'UniformOutput', 0);
% 
%     PCell = cellfun(@(x) [x' analyVar.weightPeak], PCell, 'UniformOutput', 0);
%     
%     %%% Use ceil(j/avgDataset.numTweezers) since each tweezer within the
%     %%% same image should correlate to the same independent variable
%     safeVal = regexprep(num2str(avgDataset.imagevcoAtom(k),'%.12g'),'[^a-zA-Z0-9_\-\.]','_');
% 
%     paramFile = fullfile(avgDataset.avgDir, ...
%         ['AvgFit_' analyVar.avgScanParamField '_' safeVal '_' ...
%         analyVar.sampleType analyVar.fitModel analyVar.paramFitFilename ...
%         analyVar.paramFitFileExt]);
% 
%     dlmwrite(paramFile, PCell, '\t');
% 
%     avgDataset.paramFitFiles{j} = paramFile;
% 
% end
% 
% save(fullfile(avgDataset.avgDir,'avgDataset.mat'),'avgDataset');
% 
% fclose('all');
% 
% fprintf('Averaged cloud fitting completed.\n\n');
% 
% end



function avgDataset = imagefit_NumDistFit_Averaged(varargin)
% Fit averaged OD images.
%
% Modes:
%
%   analyVar.AveragedFitMode = 'scanParameter'
%       Existing behavior:
%       fit each averaged scan-parameter/tweezer image.
%
%   analyVar.AveragedFitMode = 'allImagesByTweezer'
%       New behavior:
%       collect all individual saved OD images for each tweezer,
%       average them, and fit one image per tweezer.

%% Load variables and file data
if nargin == 0

    analyVar = AnalysisVariables;

    indivDataset = get_indiv_batch_data(analyVar);

    load( ...
        fullfile( ...
            analyVar.analyOutDir, ...
            analyVar.avgOutSubDir, ...
            'avgDataset.mat'), ...
        'avgDataset');

else

    analyVar = varargin{1};
    indivDataset = varargin{2};
    avgDataset = varargin{3};

end

%% Select fitting mode
if isfield(analyVar,'AveragedFitMode') && ...
        ~isempty(analyVar.AveragedFitMode)

    fitMode = analyVar.AveragedFitMode;

else

    fitMode = 'scanParameter';

end

validModes = {
    'scanParameter'
    'allImagesByTweezer'
    };

if ~any(strcmpi(fitMode,validModes))

    error( ...
        ['Unknown analyVar.AveragedFitMode: "%s"\n' ...
         'Valid choices are:\n' ...
         '  ''scanParameter''\n' ...
         '  ''allImagesByTweezer'''], ...
        fitMode);

end

%% Make sure output directory exists
if ~isfield(avgDataset,'avgDir') || isempty(avgDataset.avgDir)

    avgDataset.avgDir = fullfile( ...
        analyVar.analyOutDir, ...
        analyVar.avgOutSubDir);

end

if ~exist(avgDataset.avgDir,'dir')
    mkdir(avgDataset.avgDir);
end

%% Run selected fitting path
switch lower(fitMode)

    case lower('scanParameter')

        fprintf('\n');
        fprintf('Fitting averaged OD images by scan parameter...\n');

        avgDataset = fit_scan_parameter_images( ...
            analyVar, ...
            indivDataset, ...
            avgDataset);

    case lower('allImagesByTweezer')

        fprintf('\n');
        fprintf('Averaging all individual images by tweezer...\n');

        avgDataset = fit_all_images_by_tweezer( ...
            analyVar, ...
            indivDataset, ...
            avgDataset);

end

%% Record active mode
avgDataset.AveragedFitMode = fitMode;

%% Save updated avgDataset
save( ...
    fullfile(avgDataset.avgDir,'avgDataset.mat'), ...
    'avgDataset', ...
    '-v7.3');

fclose('all');

fprintf('\nAveraged cloud fitting completed.\n\n');

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avgDataset = fit_scan_parameter_images( ...
    analyVar, indivDataset, avgDataset)
% Existing behavior:
% fit each parameter/tweezer averaged image separately.

numFitsToRun = avgDataset.CounterAtom;

avgDataset.paramFitFiles = cell(numFitsToRun,1);
avgDataset.All_PCell = cell(numFitsToRun,1);

for j = 1:numFitsToRun

    fprintf( ...
        'Fitting averaged image %d of %d\n', ...
        j,numFitsToRun);

    %% Recover original-image provenance
    if isfield(avgDataset,'sourceInfo') && ...
            numel(avgDataset.sourceInfo) >= j && ...
            ~isempty(avgDataset.sourceInfo{j})

        src = avgDataset.sourceInfo{j};

        basenameNum = src.basenameNums(1);
        k = src.imageNums(1);

        if isfield(src,'tweezerNums') && ...
                ~isempty(src.tweezerNums)

            tweezerNum = src.tweezerNums(1);

        else

            tweezerNum = get_tweezer_number_from_index( ...
                j,avgDataset.numTweezers);

        end

    else

        basenameNum = 1;
        k = 1;

        tweezerNum = get_tweezer_number_from_index( ...
            j,avgDataset.numTweezers);

    end

    %% Load averaged image
    if isfield(avgDataset,'avgODFiles') && ...
            numel(avgDataset.avgODFiles) >= j && ...
            exist(avgDataset.avgODFiles{j},'file')

        OD_Image_Raw = readmatrix(avgDataset.avgODFiles{j});

    elseif isfield(avgDataset,'All_OD_Image') && ...
            numel(avgDataset.All_OD_Image) >= j

        OD_Image_Raw = avgDataset.All_OD_Image{j};

    elseif isfield(avgDataset,'avgODImages') && ...
            numel(avgDataset.avgODImages) >= j

        OD_Image_Raw = avgDataset.avgODImages{j};

    else

        error('Could not locate averaged OD image %d.',j);

    end

    %% Fit image
    [PCell,OD_Image_FitInput] = fit_one_average_image( ...
        analyVar, ...
        indivDataset, ...
        OD_Image_Raw, ...
        basenameNum, ...
        k);

    avgDataset.All_PCell{j} = PCell;

    if ~isfield(avgDataset,'All_OD_Image') || ...
            numel(avgDataset.All_OD_Image) < j

        avgDataset.All_OD_Image{j,1} = OD_Image_Raw;

    end

    %% Determine parameter value for filename
    parameterValue = get_scan_parameter_value( ...
        avgDataset,j,k);

    safeVal = make_safe_numeric_string(parameterValue);

    paramFile = fullfile( ...
        avgDataset.avgDir, ...
        ['AvgFit_' ...
         analyVar.avgScanParamField '_' ...
         safeVal '_' ...
         sprintf('Tweezer%03d_',tweezerNum) ...
         analyVar.sampleType ...
         analyVar.fitModel ...
         analyVar.paramFitFilename ...
         analyVar.paramFitFileExt]);

    write_pcell_file(paramFile,PCell);

    avgDataset.paramFitFiles{j} = paramFile;

    %% Preserve the exact image sent into the fit
    avgDataset.fitInputImages{j,1} = OD_Image_FitInput;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avgDataset = fit_all_images_by_tweezer( ...
    analyVar, indivDataset, avgDataset)
% Collect all individual saved OD images for each tweezer, average them,
% and fit one image per tweezer.

if ~isfield(avgDataset,'numTweezers') || ...
        isempty(avgDataset.numTweezers)

    error('avgDataset.numTweezers is missing.');

end

numTweezers = avgDataset.numTweezers;

if ~isfield(avgDataset,'sourceInfo') || ...
        isempty(avgDataset.sourceInfo)

    error([ ...
        'avgDataset.sourceInfo is required for allImagesByTweezer mode.\n' ...
        'Rebuild avgDataset with the version of ' ...
        'imagefit_BuildAveragedScans that saves sourceInfo.']);

end

%% Initialize new structure
tweezerFit = struct;

tweezerFit.mode = 'allImagesByTweezer';
tweezerFit.numTweezers = numTweezers;

tweezerFit.images = cell(numTweezers,1);
tweezerFit.fitInputImages = cell(numTweezers,1);

tweezerFit.imageFiles = cell(numTweezers,1);
tweezerFit.paramFitFiles = cell(numTweezers,1);

tweezerFit.PCell = cell(numTweezers,1);

tweezerFit.sourceODFiles = cell(numTweezers,1);
tweezerFit.sourceBasenameNums = cell(numTweezers,1);
tweezerFit.sourceImageNums = cell(numTweezers,1);

tweezerFit.numImagesAveraged = zeros(numTweezers,1);

tweezerFit.referenceBasenameNums = nan(numTweezers,1);
tweezerFit.referenceImageNums = nan(numTweezers,1);

%% Loop through tweezers
for tweezerNum = 1:numTweezers

    fprintf( ...
        'Collecting individual images for tweezer %d of %d\n', ...
        tweezerNum,numTweezers);

    [odFiles,basenameNums,imageNums] = ...
        collect_individual_od_files_for_tweezer( ...
            avgDataset,tweezerNum);

    if isempty(odFiles)

        warning( ...
            'No individual OD files were found for tweezer %d.', ...
            tweezerNum);

        continue;

    end

    %% Remove repeated file references
    %
    % The same source file should normally appear in only one
    % parameter/tweezer group, but unique prevents accidental duplication.

    [odFiles,uniqueIndex] = unique(odFiles,'stable');

    basenameNums = basenameNums(uniqueIndex);
    imageNums = imageNums(uniqueIndex);

    %% Read individual images
    imageCell = cell(numel(odFiles),1);

    validImage = false(numel(odFiles),1);

    referenceSize = [];

    for imageCounter = 1:numel(odFiles)

        thisFile = odFiles{imageCounter};

        if ~exist(thisFile,'file')

            warning('Missing OD image file:\n%s',thisFile);
            continue;

        end

        thisImage = readmatrix(thisFile);

        if isempty(referenceSize)

            referenceSize = size(thisImage);

        elseif ~isequal(size(thisImage),referenceSize)

            warning([ ...
                'Skipping OD image because its size does not match ' ...
                'the other tweezer images:\n%s'], ...
                thisFile);

            continue;

        end

        imageCell{imageCounter} = double(thisImage);
        validImage(imageCounter) = true;

    end

    %% Keep only valid images
    imageCell = imageCell(validImage);
    odFiles = odFiles(validImage);
    basenameNums = basenameNums(validImage);
    imageNums = imageNums(validImage);

    if isempty(imageCell)

        warning( ...
            'No readable OD images remained for tweezer %d.', ...
            tweezerNum);

        continue;

    end

    %% Build image stack
    imageSize = size(imageCell{1});

    imageStack = nan( ...
        imageSize(1), ...
        imageSize(2), ...
        numel(imageCell));

    for imageCounter = 1:numel(imageCell)

        imageStack(:,:,imageCounter) = ...
            imageCell{imageCounter};

    end

    %% Average every individual image for this tweezer
    averageImage = mean(imageStack,3,'omitnan');

    %% Save averaged image
    averageImageFile = fullfile( ...
        avgDataset.avgDir, ...
        sprintf( ...
            'AllImages_Average_Tweezer%03d.txt', ...
            tweezerNum));

    writematrix( ...
        averageImage, ...
        averageImageFile, ...
        'Delimiter','tab');

    %% Use the first source image as fitting-geometry reference
    basenameNum = basenameNums(1);
    k = imageNums(1);

    %% Fit the all-images average
    fprintf( ...
        'Fitting all-images average for tweezer %d using %d images\n', ...
        tweezerNum,numel(imageCell));

    [PCell,OD_Image_FitInput] = fit_one_average_image( ...
        analyVar, ...
        indivDataset, ...
        averageImage, ...
        basenameNum, ...
        k);

    %% Save fit coefficients
    paramFile = fullfile( ...
        avgDataset.avgDir, ...
        ['AllImagesFit_' ...
         sprintf('Tweezer%03d_',tweezerNum) ...
         analyVar.sampleType ...
         analyVar.fitModel ...
         analyVar.paramFitFilename ...
         analyVar.paramFitFileExt]);

    write_pcell_file(paramFile,PCell);

    %% Store results
    tweezerFit.images{tweezerNum} = averageImage;
    tweezerFit.fitInputImages{tweezerNum} = OD_Image_FitInput;

    tweezerFit.imageFiles{tweezerNum} = averageImageFile;
    tweezerFit.paramFitFiles{tweezerNum} = paramFile;

    tweezerFit.PCell{tweezerNum} = PCell;

    tweezerFit.sourceODFiles{tweezerNum} = odFiles;
    tweezerFit.sourceBasenameNums{tweezerNum} = basenameNums;
    tweezerFit.sourceImageNums{tweezerNum} = imageNums;

    tweezerFit.numImagesAveraged(tweezerNum) = numel(imageCell);

    tweezerFit.referenceBasenameNums(tweezerNum) = basenameNum;
    tweezerFit.referenceImageNums(tweezerNum) = k;

end

%% Save new fitting collection
avgDataset.tweezerAllImageFit = tweezerFit;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [odFiles,basenameNums,imageNums] = ...
    collect_individual_od_files_for_tweezer( ...
        avgDataset,requestedTweezer)
% Recover every individual saved OD file belonging to one tweezer.

odFiles = {};
basenameNums = [];
imageNums = [];

for j = 1:numel(avgDataset.sourceInfo)

    src = avgDataset.sourceInfo{j};

    if isempty(src) || ...
            ~isfield(src,'tweezerNums') || ...
            ~isfield(src,'odFiles')

        continue;

    end

    sourceTweezers = src.tweezerNums(:);

    sourceODFiles = src.odFiles;

    if ischar(sourceODFiles) || isstring(sourceODFiles)
        sourceODFiles = cellstr(sourceODFiles);
    end

    sourceODFiles = sourceODFiles(:);

    if isfield(src,'basenameNums')
        sourceBasenames = src.basenameNums(:);
    else
        sourceBasenames = nan(numel(sourceODFiles),1);
    end

    if isfield(src,'imageNums')
        sourceImages = src.imageNums(:);
    else
        sourceImages = nan(numel(sourceODFiles),1);
    end

    numSourceRows = min([
        numel(sourceTweezers)
        numel(sourceODFiles)
        numel(sourceBasenames)
        numel(sourceImages)
        ]);

    sourceTweezers = sourceTweezers(1:numSourceRows);
    sourceODFiles = sourceODFiles(1:numSourceRows);
    sourceBasenames = sourceBasenames(1:numSourceRows);
    sourceImages = sourceImages(1:numSourceRows);

    matchingRows = sourceTweezers == requestedTweezer;

    odFiles = [
        odFiles
        sourceODFiles(matchingRows)
        ];

    basenameNums = [
        basenameNums
        sourceBasenames(matchingRows)
        ];

    imageNums = [
        imageNums
        sourceImages(matchingRows)
        ];

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PCell,OD_Image_FitInput] = fit_one_average_image( ...
    analyVar, ...
    indivDataset, ...
    OD_Image_Raw, ...
    basenameNum, ...
    k)
% Apply the existing fitting procedure to one averaged image.

OD_Image_FitInput = double(OD_Image_Raw);

%% Optional smoothing
if analyVar.fitSmoothOD

    OD_Image_FitInput = analyVar.smoothFilt( ...
        OD_Image_FitInput, ...
        analyVar.smoothFiltMat);

end

%% Construct fit windows
if isfield(analyVar,'UseTweezer') && ...
        analyVar.UseTweezer == 1

    % The saved tweezer image is already the entire fit window.
    roiWin_Index = {
        true(size(OD_Image_FitInput))
        };

    if size(OD_Image_FitInput,1) ~= ...
            size(OD_Image_FitInput,2)

        error([ ...
            'The current fitting code expects a square tweezer image, ' ...
            'but received size %d x %d.'], ...
            size(OD_Image_FitInput,1), ...
            size(OD_Image_FitInput,2));

    end

    fitWinSize = size(OD_Image_FitInput,1);

    OD_Fit_ImageCell = {
        OD_Image_FitInput
        };

else

    if isfield(indivDataset{basenameNum},'roiWin_Index')

        roiWin_Index = ...
            indivDataset{basenameNum}.roiWin_Index;

    else

        error( ...
            'indivDataset{%d}.roiWin_Index is missing.', ...
            basenameNum);

    end

    fitWinSize = analyVar.funcFitWin(basenameNum);

    fitLogicCell = cellfun( ...
        analyVar.fitWinLogicInd, ...
        roiWin_Index, ...
        'UniformOutput',false);

    OD_Fit_ImageCell = cellfun( ...
        @(fitLogic) reshape( ...
            OD_Image_FitInput(fitLogic), ...
            [fitWinSize fitWinSize]), ...
        fitLogicCell, ...
        'UniformOutput',false);

end

%% Initial guesses
InitGuess = get_fit_params( ...
    analyVar, ...
    OD_Image_FitInput, ...
    OD_Fit_ImageCell, ...
    indivDataset, ...
    basenameNum, ...
    k);

%% Bounds
[lowBndVec,upBndVec] = ...
    get_fit_bounds(analyVar,basenameNum);

nFits = numel(InitGuess);

lowBndCell = cell(size(InitGuess));
upBndCell = cell(size(InitGuess));

for q = 1:nFits

    lowBndCell{q} = lowBndVec(:);
    upBndCell{q} = upBndVec(:);

end

%% Weighting
errCell = cellfun( ...
    @(imageData) get_OD_weight( ...
        analyVar.weightPeak,imageData), ...
    OD_Fit_ImageCell, ...
    'UniformOutput',false);

%% Fit coordinate grid
[Xgrid,Ygrid] = meshgrid(1:fitWinSize);

optimOpt = optimset( ...
    'Display','off', ...
    'FinDiffType','central', ...
    'TolFun',1e-9, ...
    'TolX',1e-9);

%% Fit
PCell = cellfun( ...
    @(fitImage,fitError,initialGuess,lowerBounds,upperBounds) ...
        lsqcurvefit( ...
            str2func(analyVar.fitModel), ...
            initialGuess, ...
            [Xgrid(:),Ygrid(:),fitError(:)], ...
            fitImage(:)./fitError(:), ...
            lowerBounds, ...
            upperBounds, ...
            optimOpt), ...
    OD_Fit_ImageCell, ...
    errCell, ...
    InitGuess, ...
    lowBndCell, ...
    upBndCell, ...
    'UniformOutput',false);

%% Append the weighting parameter
PCell = cellfun( ...
    @(fitCoefficients) [ ...
        fitCoefficients(:)' ...
        analyVar.weightPeak], ...
    PCell, ...
    'UniformOutput',false);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_pcell_file(paramFile,PCell)
% Write one row per fit region.

if isempty(PCell)

    warning('No fit coefficients were generated for:\n%s',paramFile);
    return;

end

PCellRows = cellfun( ...
    @(row) row(:)', ...
    PCell, ...
    'UniformOutput',false);

PMatrix = vertcat(PCellRows{:});

writematrix( ...
    PMatrix, ...
    paramFile, ...
    'Delimiter','tab');

end



function tweezerNum = get_tweezer_number_from_index( ...
    imageIndex,numTweezers)

tweezerNum = mod(imageIndex-1,numTweezers)+1;

end


function parameterValue = get_scan_parameter_value( ...
    avgDataset,j,k)

if isfield(avgDataset,'paramVals') && ...
        numel(avgDataset.paramVals) >= j

    parameterValue = avgDataset.paramVals(j);

elseif isfield(avgDataset,'imagevcoAtom')

    if numel(avgDataset.imagevcoAtom) >= j

        parameterValue = avgDataset.imagevcoAtom(j);

    elseif numel(avgDataset.imagevcoAtom) >= k

        parameterValue = avgDataset.imagevcoAtom(k);

    else

        parameterValue = j;

    end

else

    parameterValue = j;

end

end

function safeVal = make_safe_numeric_string(value)

safeVal = regexprep( ...
    num2str(value,'%.12g'), ...
    '[^a-zA-Z0-9_\-\.]', ...
    '_');

end