function avgDataset = imagefit_NumDistFit_AveragedV2(varargin)
% Fit averaged OD images for each scan-parameter/tweezer combination.
%
% Workflow:
%   1. Load analysis variables and avgDataset from imagefit_BuildAveragedScans.
%   2. Locate each averaged OD image.
%   3. Fit the averaged image.
%   4. Save fit parameters and the images used for fitting.
%   5. Save the updated avgDataset.
%
% This function uses the scan-parameter fitting workflow only.
%
% Main outputs added to avgDataset:
%   All_OD_Image       - averaged OD images used for fitting
%   All_PCell          - fitted parameter cells
%   fitInputImages     - exact images passed into the fitting routine
%   paramFitFiles      - files containing fitted parameters

%% ================================================================
%  LOAD VARIABLES AND DATA
%  ================================================================

if nargin == 0

    % Normal standalone use
    analyVar = AnalysisVariables;

    indivDataset = get_indiv_batch_data(analyVar);

    load( ...
        fullfile( ...
            analyVar.analyOutDir, ...
            analyVar.avgOutSubDir, ...
            'avgDataset.mat'), ...
        'avgDataset');

else

    % Allow existing variables to be passed in
    analyVar = varargin{1};
    indivDataset = varargin{2};
    avgDataset = varargin{3};

end


%% ================================================================
%  GRAB AVGDATASET OUTPUT DIRECTORY
%  ================================================================


avgDataset.avgDir = fullfile( ...
        analyVar.analyOutDir, ...
        analyVar.avgOutSubDir);


%% ================================================================
%  FIT ALL AVERAGED SCAN-PARAMETER IMAGES
%  ================================================================

fprintf('\n');
fprintf('Fitting averaged OD images by scan parameter...\n\n');

avgDataset = fit_scan_parameter_images( ...
    analyVar, ...
    indivDataset, ...
    avgDataset);


%% ================================================================
%  SAVE UPDATED DATASET
%  ================================================================

save( ...
    fullfile(avgDataset.avgDir,'avgDataset.mat'), ...
    'avgDataset', ...
    '-v7.3');

fclose('all');

fprintf('\nAveraged cloud fitting completed.\n\n');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avgDataset = fit_scan_parameter_images( ...
    analyVar,indivDataset,avgDataset)
% Fit each averaged scan-parameter/tweezer image.
%
% For each averaged image this function:
%   - Finds the source image information.
%   - Loads the averaged OD image.
%   - Performs the 2D fit.
%   - Saves the fitted parameters.
%   - Stores the image and fit results in avgDataset.

%% ================================================================
%  INITIALIZE OUTPUT ARRAYS
%  ================================================================

numFitsToRun = avgDataset.CounterAtom;

avgDataset.paramFitFiles = cell(numFitsToRun,1);
avgDataset.All_PCell = cell(numFitsToRun,1);
avgDataset.All_OD_Image = cell(numFitsToRun,1);
avgDataset.fitInputImages = cell(numFitsToRun,1);


%% ================================================================
%  FIT EACH AVERAGED IMAGE
%  ================================================================

for j = 1:numFitsToRun

    fprintf( ...
        'Fitting averaged image %d of %d\n', ...
        j,numFitsToRun);


    %% ------------------------------------------------------------
    %  FIND SOURCE IMAGE INFORMATION
    %  ------------------------------------------------------------

    if isfield(avgDataset,'sourceInfo') && ...
            numel(avgDataset.sourceInfo) >= j && ...
            ~isempty(avgDataset.sourceInfo{j})

        src = avgDataset.sourceInfo{j};

        basenameNum = src.basenameNums(1);
        k = src.imageNums(1);

        % Use the stored tweezer number when available.
        if isfield(src,'tweezerNums') && ...
                ~isempty(src.tweezerNums)

            tweezerNum = src.tweezerNums(1);

        else

            tweezerNum = get_tweezer_number_from_index( ...
                j, ...
                avgDataset.numTweezers);

        end

    else

        % Fallback if source information was not saved.
        basenameNum = 1;
        k = 1;

        tweezerNum = get_tweezer_number_from_index( ...
            j, ...
            avgDataset.numTweezers);

    end


    %% ------------------------------------------------------------
    %  LOAD THE AVERAGED OD IMAGE
    %  ------------------------------------------------------------
    %
    % Preferred source:
    %   IndivTwzrAvgODFiles
    %
    % Fallbacks are kept for compatibility with older avgDatasets.

    if isfield(avgDataset,'IndivTwzrAvgODFiles') && ...
            numel(avgDataset.IndivTwzrAvgODFiles) >= j && ...
            exist(avgDataset.IndivTwzrAvgODFiles{j},'file')

        OD_Image_Raw = readmatrix( ...
            avgDataset.IndivTwzrAvgODFiles{j});

    elseif isfield(avgDataset,'All_OD_Image') && ...
            numel(avgDataset.All_OD_Image) >= j && ...
            ~isempty(avgDataset.All_OD_Image{j})

        OD_Image_Raw = avgDataset.All_OD_Image{j};

    elseif isfield(avgDataset,'IndivTwzrAvgODImages') && ...
            numel(avgDataset.IndivTwzrAvgODImages) >= j && ...
            ~isempty(avgDataset.IndivTwzrAvgODImages{j})

        OD_Image_Raw = avgDataset.IndivTwzrAvgODImages{j};

    else

        error( ...
            'Could not locate averaged OD image %d.', ...
            j);

    end


    % Make sure the image is actually usable before fitting.
    if isempty(OD_Image_Raw)

        error( ...
            'Averaged OD image %d is empty.', ...
            j);

    end

    fprintf( ...
        '  Loaded averaged image %d: size = [%d %d], empty = %d\n', ...
        j, ...
        size(OD_Image_Raw,1), ...
        size(OD_Image_Raw,2), ...
        isempty(OD_Image_Raw));


    %% ------------------------------------------------------------
    %  FIT THE AVERAGED IMAGE
    %  ------------------------------------------------------------

    [PCell,OD_Image_FitInput] = fit_one_average_image( ...
        analyVar, ...
        indivDataset, ...
        OD_Image_Raw, ...
        basenameNum, ...
        k);


    %% ------------------------------------------------------------
    %  STORE FIT RESULTS
    %  ------------------------------------------------------------

    % Fitted parameter cell.
    avgDataset.All_PCell{j,1} = PCell;

    % Exact averaged image used for the fit.
    avgDataset.All_OD_Image{j,1} = OD_Image_Raw;

    % Exact image passed into the fitting routine.
    % This may differ from All_OD_Image if smoothing is enabled.
    avgDataset.fitInputImages{j,1} = OD_Image_FitInput;


    %% ------------------------------------------------------------
    %  CREATE FIT PARAMETER FILENAME
    %  ------------------------------------------------------------

    parameterValue = get_scan_parameter_value( ...
        avgDataset, ...
        j, ...
        k);

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


    %% ------------------------------------------------------------
    %  SAVE FIT PARAMETERS
    %  ------------------------------------------------------------

    write_pcell_file( ...
        paramFile, ...
        PCell);

    avgDataset.paramFitFiles{j} = paramFile;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PCell,OD_Image_FitInput] = fit_one_average_image( ...
    analyVar, ...
    indivDataset, ...
    OD_Image_Raw, ...
    basenameNum, ...
    k)
% Perform the actual 2D fit on one averaged OD image.
%
% This function handles:
%   1. Optional smoothing
%   2. Fit-window construction
%   3. Initial guesses
%   4. Fit bounds
%   5. Weighting
%   6. Nonlinear least-squares fitting
%   7. Adding the weighting parameter to PCell


%% ================================================================
%  PREPARE FIT INPUT IMAGE
%  ================================================================

OD_Image_FitInput = double(OD_Image_Raw);

% Optional smoothing used as the fitting input.
if analyVar.fitSmoothOD

    OD_Image_FitInput = analyVar.smoothFilt( ...
        OD_Image_FitInput, ...
        analyVar.smoothFiltMat);

end


%% ================================================================
%  CONSTRUCT FIT WINDOWS
%  ================================================================

if isfield(analyVar,'UseTweezer') && ...
        analyVar.UseTweezer == 1

    % --------------------------------------------------------------
    % Tweezer mode:
    % Each saved averaged image is already one complete fit window.
    % --------------------------------------------------------------

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

    % --------------------------------------------------------------
    % Full-cloud mode:
    % Use the saved ROI windows from indivDataset.
    % --------------------------------------------------------------

    if isfield(indivDataset{basenameNum},'roiWin_Index')

        roiWin_Index = ...
            indivDataset{basenameNum}.roiWin_Index;

    else

        error( ...
            'indivDataset{%d}.roiWin_Index is missing.', ...
            basenameNum);

    end

    fitWinSize = analyVar.funcFitWin(basenameNum);

    % Convert each ROI definition into a logical fit mask.
    fitLogicCell = cellfun( ...
        analyVar.fitWinLogicInd, ...
        roiWin_Index, ...
        'UniformOutput',false);

    % Extract the corresponding pixels from the averaged image.
    OD_Fit_ImageCell = cellfun( ...
        @(fitLogic) reshape( ...
            OD_Image_FitInput(fitLogic), ...
            [fitWinSize fitWinSize]), ...
        fitLogicCell, ...
        'UniformOutput',false);

end


%% ================================================================
%  INITIAL FIT GUESSES
%  ================================================================

InitGuess = get_fit_params( ...
    analyVar, ...
    OD_Image_FitInput, ...
    OD_Fit_ImageCell, ...
    indivDataset, ...
    basenameNum, ...
    k);


%% ================================================================
%  FIT PARAMETER BOUNDS
%  ================================================================

[lowBndVec,upBndVec] = ...
    get_fit_bounds( ...
        analyVar, ...
        basenameNum);

nFits = numel(InitGuess);

lowBndCell = cell(size(InitGuess));
upBndCell = cell(size(InitGuess));

for q = 1:nFits

    lowBndCell{q} = lowBndVec(:);
    upBndCell{q} = upBndVec(:);

end


%% ================================================================
%  FITTING WEIGHT
%  ================================================================

errCell = cellfun( ...
    @(imageData) ...
        get_OD_weight( ...
            analyVar.weightPeak, ...
            imageData), ...
    OD_Fit_ImageCell, ...
    'UniformOutput',false);


%% ================================================================
%  FIT COORDINATE GRID
%  ================================================================

[Xgrid,Ygrid] = meshgrid(1:fitWinSize);


%% ================================================================
%  OPTIMIZATION SETTINGS
%  ================================================================

optimOpt = optimset( ...
    'Display','off', ...
    'FinDiffType','central', ...
    'TolFun',1e-9, ...
    'TolX',1e-9);

        %% ================================================================
        %  PERFORM 1D FIT
        %  ================================================================
        
 if isfield(analyVar,'dimenReduc1D') && ...
                analyVar.dimenReduc1D == 1
        
            PCell = fit_integrated_1D( ...
                analyVar, ...
                OD_Fit_ImageCell, ...
                InitGuess, ...
                lowBndVec, ...
                upBndVec, ...
                optimOpt);
        
else
        %% ================================================================
        %  PERFORM NONLINEAR FIT
        %  ================================================================
        
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
end


%% ================================================================
%  APPEND FIT WEIGHT TO PARAMETER CELL
%  ================================================================

PCell = cellfun( ...
    @(fitCoefficients) ...
        [ ...
        fitCoefficients(:)' ...
        analyVar.weightPeak], ...
    PCell, ...
    'UniformOutput',false);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_pcell_file(paramFile,PCell)
% Write fitted parameters to a tab-delimited text file.
%
% Each row corresponds to one fitted region.
% The final column contains the weighting parameter.


if isempty(PCell)

    warning( ...
        'No fit coefficients were generated for:\n%s', ...
        paramFile);

    return;

end


% Convert each parameter row into a row vector.
PCellRows = cellfun( ...
    @(row) row(:)', ...
    PCell, ...
    'UniformOutput',false);


% Combine all fit rows into one numeric matrix.
PMatrix = vertcat(PCellRows{:});


% Save as tab-delimited text.
writematrix( ...
    PMatrix, ...
    paramFile, ...
    'Delimiter','tab');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tweezerNum = get_tweezer_number_from_index( ...
    imageIndex,numTweezers)
% Convert the sequential averaged-image index into a tweezer number.
%
% Example for 3 tweezers:
%   image 1 -> tweezer 1
%   image 2 -> tweezer 2
%   image 3 -> tweezer 3
%   image 4 -> tweezer 1
%   image 5 -> tweezer 2
%   ...


tweezerNum = ...
    mod(imageIndex-1,numTweezers)+1;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameterValue = get_scan_parameter_value( ...
    avgDataset,j,k)
% Get the scan-parameter value associated with averaged image j.
%
% Preferred field:
%   avgDataset.paramVals
%
% Fallback:
%   avgDataset.imagevcoAtom


if isfield(avgDataset,'paramVals') && ...
        numel(avgDataset.paramVals) >= j

    parameterValue = ...
        avgDataset.paramVals(j);

elseif isfield(avgDataset,'imagevcoAtom')

    if numel(avgDataset.imagevcoAtom) >= j

        parameterValue = ...
            avgDataset.imagevcoAtom(j);

    elseif numel(avgDataset.imagevcoAtom) >= k

        parameterValue = ...
            avgDataset.imagevcoAtom(k);

    else

        parameterValue = j;

    end

else

    parameterValue = j;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function safeVal = make_safe_numeric_string(value)
% Convert a numeric scan parameter into a filename-safe string.


safeVal = regexprep( ...
    num2str(value,'%.12g'), ...
    '[^a-zA-Z0-9_\-\.]', ...
    '_');

end