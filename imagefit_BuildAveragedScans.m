function avgDataset = imagefit_BuildAveragedScans(varargin)

%% Load variables and file data
if nargin == 0
    analyVar = AnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
end

avgDir = fullfile(analyVar.analyOutDir, analyVar.avgOutSubDir);
if ~exist(avgDir,'dir')
    mkdir(avgDir);
end

fprintf('\nBuilding averaged OD images across matching scans...\n');

% Load tweezer ROIs

load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
numTweezers = size(tweezerROI.centersXY,1);

% From each individual file in indivDatset, Collect all values into a single list
% to scan over later when creating the properly averaged lists:

allParamVals    = [];   %   scanned parameter values (Or Dummy Scan file ID)
allTweezerNums  = [];   %   tweezer number this entry comes from
allODImages       = {};   %   saved image in the OD file (subracted counts from bkg image)
allRawAtomImgFiles      = {};   %   collection of raw atom image files
allODFiles      = {};   %   name of image file used
allTotODCounts  = [];   %   number of counts in tweezer ROI (again subtracted)
allTotCountsImg1Raw = [];     % number of raw counts in the ROI 

allBasenameNums = [];   %   collecting the basename number for each entry
allImageNums    = [];   %   collecting the appearance number of the entry within the file

scanIDs = analyVar.meanListVar;
allScanIDs = [];

for basenameNum = 1:analyVar.numBasenamesAtom

    paramVals = get_average_scan_param_values(analyVar, indivDataset, basenameNum); % will store the independent variable to plot against based scanned param or given ID

    for k = 1:indivDataset{basenameNum}.CounterAtom

        for tweezerNum = 1:numTweezers
            
            % grab OD image
            odFile = [analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                sprintf('_Tweezer%03d',tweezerNum) ...
                analyVar.ODimageFilename];

            if ~exist(odFile,'file')
                error('Missing OD image file:\n%s', odFile);
            end
    
            OD = dlmread(odFile);


            allParamVals(end+1,1)    = paramVals(k);    %The corresponding imagevcoAtom value of entry
            allTweezerNums(end+1,1)  = tweezerNum;  %Tweezer number entry
            allODImages{end+1,1}       = OD;        %OD image matrix
            allRawAtomImgFiles{end+1,1}      = char(indivDataset{basenameNum}.fileAtom(k));     %Raw atom image file name
            allODFiles{end+1,1}      = odFile;      %OD image file
            allBasenameNums(end+1,1) = basenameNum;    %basename number of entry
            allImageNums(end+1,1)    = k;   %entry number within the basename num (file)
            allTotODCounts(end+1,1) = indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum);                   % PCA Bkg Subtracted Counts (May be normalized if toggled)
            allTotCountsImg1Raw(end+1,1) = indivDataset{basenameNum}.TotalCountsImg1Raw(k,tweezerNum);     % Raw Counts in atom image
            allScanIDs(end+1,1) = scanIDs(basenameNum);
        end
    end
end

%    %% Group by parameter value
%    Will group by common independent values (avg across multiple files) or
%    group by all scans in a single file (avg a single dummy scan)
%    uniqueVals = unique(allParamVals);
%    
%    %% Initialize avgDataset
%    
%    % Any time there is a reference like: indivDataset{j}.whatever
%    Where the source of the data in avg dataset is needed it will be
%    stored in the source struct within avgdataset
%    %        src = avgDataset.sourceInfo{j};
%    %    
%    %        basenameNum = src.basenameNums(1);
%    %        k           = src.imageNums(1);
%    %        tweezerNum  = src.tweezerNums(1);
%    %        
%    %        value = indivDataset{basenameNum}.whatever(k);


%% Group by scan parameter value within AvgDataset

allParamVals   = allParamVals(:);
allTweezerNums = allTweezerNums(:);
allTotODCounts = allTotODCounts(:);
allTotImg1Counts = allTotCountsImg1Raw(:);

%% Find unique paramVals

uniqueVals = unique(allParamVals,'sorted');
uniqueScanIDs = unique(allScanIDs,'sorted');


%% Initialize avgDataset Individual entry for each tweezer in a file

avgDataset = struct;

avgDataset.avgDir       = avgDir;
avgDataset.numTweezers  = numTweezers;
avgDataset.TweezerROI   = tweezerROI;

avgDataset.CounterAtom  = 0;
avgDataset.tweezerNums  = [];
avgDataset.imagevcoAtom = [];

avgDataset.IndivTwzrAvgODFiles   = {};
avgDataset.IndivTwzrAvgODImages  = {};
avgDataset.sourceFiles  = {};
avgDataset.sourceInfo   = {};

avgDataset.IndivTwzrRawCounts       = {};
avgDataset.IndivTwzrTotODCounts       = {};

avgDataset.IndivTwzrMeanRawODCounts   = [];
avgDataset.IndivTwzrStdRawODCounts    = [];

avgDataset.IndivTwzrMeanODCounts      = [];
avgDataset.IndivTwzrStdODCounts       = [];

avgDataset.IndivTwzrNumAveraged = [];

% Compatibility fields
avgDataset.roiWin_Index = indivDataset{1}.roiWin_Index;

% Record grouping information
avgDataset.uniqueImagevcoAtom  = uniqueVals;

%% Create Collection of the mean values of individual tweezers in same group


%% Build averaged parameter/tweezer entries

avgDatasetCounter = 0;

for p = 1:numel(uniqueVals)

    groupValue = uniqueVals(p);

    for s = 1:numel(uniqueScanIDs)

        groupScanID = uniqueScanIDs(s);

        for tweezerNum = 1:numTweezers

            idx = find( ...
                allParamVals == groupValue & ...
                allScanIDs == groupScanID & ...
                allTweezerNums == tweezerNum);

            % This paramVal does not exist in this scan
            if isempty(idx)
                continue;
            end

        %selectedParamVals = allParamVals(idx);

        %% Grab image dimensions using first image and create image stack

        firstImage = allODImages{idx(1)};

        imageRows = size(firstImage,1);
        imageCols = size(firstImage,2);

        imgStack = nan( ...
            imageRows, ...
            imageCols, ...
            numel(idx));

        %% Sum up images of same parameters

        for n = 1:numel(idx)

            recordIndex = idx(n);
            thisImage = allODImages{recordIndex};

            if ~isequal(size(thisImage),[imageRows imageCols])
                error(['Image-size mismatch for record %d. Expected ' ...
                       '%d-by-%d but found %d-by-%d.'], ...
                    recordIndex, ...
                    imageRows,imageCols, ...
                    size(thisImage,1),size(thisImage,2));
            end

            imgStack(:,:,n) = thisImage;

        end

        %% Average images

        avgOD = mean(imgStack,3,'omitnan');

        %% Collect externally calculated counts

        totODCounts = allTotODCounts(idx);

        rawCounts = allTotImg1Counts(idx);

        %% Store averaged entry sorting information
        % One entry for every averaged parameter/tweezer image

        avgDatasetCounter = avgDatasetCounter + 1;

        avgDataset.CounterAtom = avgDatasetCounter;

        avgDataset.tweezerNums(avgDatasetCounter,1) = tweezerNum;

        avgDataset.imagevcoAtom(avgDatasetCounter,1) = groupValue;

        avgDataset.scanID(avgDatasetCounter,1) = groupScanID;

        avgDataset.IndivTwzrAvgODImages{avgDatasetCounter,1} = avgOD;

        avgDataset.sourceFiles{avgDatasetCounter,1} = allRawAtomImgFiles(idx);

        avgDataset.IndivTwzrNumAveraged(avgDatasetCounter,1) = numel(idx);

        %% Store count information

        avgDataset.IndivTwzrRawCounts{avgDatasetCounter,1} = rawCounts;

        avgDataset.IndivTwzrTotODCounts{avgDatasetCounter,1} = totODCounts;

        avgDataset.IndivTwzrMeanRawODCounts(avgDatasetCounter,1) = mean(rawCounts,'omitnan');

        avgDataset.IndivTwzrStdRawODCounts(avgDatasetCounter,1) = std(rawCounts,0,'omitnan');

        avgDataset.IndivTwzrMeanODCounts(avgDatasetCounter,1) = mean(totODCounts,'omitnan');

        avgDataset.IndivTwzrStdODCounts(avgDatasetCounter,1) = std(totODCounts,0,'omitnan');

        %% Store indivdataset information for each entry within the average

        sourceInfo = struct;

        sourceInfo.recordIndices = idx;

        sourceInfo.basenameNums = allBasenameNums(idx);

        sourceInfo.imageNums = allImageNums(idx);

        sourceInfo.tweezerNums = allTweezerNums(idx);

        sourceInfo.imagevcoAtom = groupValue;

        sourceInfo.sourceFiles = allRawAtomImgFiles(idx);

        sourceInfo.odFiles = allODFiles(idx);

        sourceInfo.rawODCounts = rawCounts;

        sourceInfo.totODCounts = totODCounts;
        
        sourceInfo.scanID = groupScanID;

        avgDataset.sourceInfo{avgDatasetCounter,1} = sourceInfo;

        %% Save averaged image

        safeVal = regexprep( ...
            num2str(groupValue,'%.12g'), ...
            '[^a-zA-Z0-9_\-\.]', ...
            '_');

        avgFile = fullfile( ...
            avgDataset.avgDir, ...
            sprintf( ...
                'AvgOD_%s_%s_Tweezer%03d.txt', ...
                analyVar.avgScanParamField, ...
                safeVal, ...
                tweezerNum));

        dlmwrite(avgFile,avgOD,'\t');

        avgDataset.IndivTwzrAvgODFiles{avgDatasetCounter,1} = ...
            avgFile;

        fprintf(['Avg entry %d: parameter group %.12g, ' ...
                 ' tweezer %d, using %d images\n'], ...
            avgDatasetCounter, ...
            groupValue, ...
            tweezerNum, ...
            numel(idx));
        end
    end
end


%% Average counts across tweezers for each scan ID and parameter value (imagevcoAtom)
%
% Each cell corresponds to one scan ID.
%
% For each scan:
%   RawCountsX{s}   = parameter values (x)
%   RawCountsY{s}   = mean raw counts (y)
%   RawCountsErr{s} = SEM of raw counts
%
%   ODCountsX{s}    = parameter values (x)
%   ODCountsY{s}    = mean OD counts (y)
%   ODCountsErr{s}  = SEM of OD counts
%
% The averaging is performed across all tweezer entries having the
% same scan ID and parameter value.

numScans = numel(uniqueScanIDs);

%% Initialize 1D cell arrays
%
% One cell per scan ID.

avgDataset.RawCountsX   = cell(numScans,1);
avgDataset.RawCountsY   = cell(numScans,1);
avgDataset.RawCountsErr = cell(numScans,1);

avgDataset.ODCountsX    = cell(numScans,1);
avgDataset.ODCountsY    = cell(numScans,1);
avgDataset.ODCountsErr  = cell(numScans,1);

avgDataset.uniqueScanIDs = uniqueScanIDs;

%% Loop over each scan ID

for s = 1:numScans

    groupScanID = uniqueScanIDs(s);

    % These will contain only parameter values that actually exist
    % in this scan.
    scanRawX   = [];
    scanRawY   = [];
    scanRawErr = [];

    scanODX    = [];
    scanODY    = [];
    scanODErr  = [];

    %% Loop over parameter values

    for p = 1:numel(uniqueVals)

        groupValue = uniqueVals(p);

        %% Find all averaged tweezer entries for this
        % parameter value and scan ID

        idx = find( ...
            avgDataset.imagevcoAtom == groupValue & ...
            avgDataset.scanID == groupScanID);

        if isempty(idx)
            continue;
        end

        %% Get the mean count for each individual tweezer
        %
        % Each element here represents one tweezer, where the
        % repeated images for that tweezer have already been
        % averaged above.

        rawMeans = avgDataset.IndivTwzrMeanRawODCounts(idx);
        odMeans  = avgDataset.IndivTwzrMeanODCounts(idx);

        %% Remove NaN/Inf values

        validRaw = isfinite(rawMeans);
        validOD  = isfinite(odMeans);

        rawMeans = rawMeans(validRaw);
        odMeans  = odMeans(validOD);

        %% Raw counts

        if ~isempty(rawMeans)

            rawMean = mean(rawMeans);

            % Standard error of the mean across tweezers
            if numel(rawMeans) > 1
                rawErr = std(rawMeans,0) / sqrt(numel(rawMeans));
            else
                rawErr = NaN;
            end

            scanRawX(end+1,1)   = groupValue;
            scanRawY(end+1,1)   = rawMean;
            scanRawErr(end+1,1) = rawErr;

        end

        %% OD counts

        if ~isempty(odMeans)

            odMean = mean(odMeans);

            % Standard error of the mean across tweezers
            if numel(odMeans) > 1
                odErr = std(odMeans,0) / sqrt(numel(odMeans));
            else
                odErr = NaN;
            end

            scanODX(end+1,1)   = groupValue;
            scanODY(end+1,1)   = odMean;
            scanODErr(end+1,1) = odErr;

        end

    end

    %% Store this scan's 1D data in the corresponding cell

    avgDataset.RawCountsX{s}   = scanRawX;
    avgDataset.RawCountsY{s}   = scanRawY;
    avgDataset.RawCountsErr{s} = scanRawErr;

    avgDataset.ODCountsX{s}    = scanODX;
    avgDataset.ODCountsY{s}    = scanODY;
    avgDataset.ODCountsErr{s}  = scanODErr;

end
%  % 
%  % % Number of groups
%  % n = numTweezers;
%  % numGroups = avgDatasetCounter / n;
%  % 
%  % % Reshape into groups
%  % MeanRawODCountGroups = reshape(avgDataset.IndivTwzrMeanRawODCounts, n, numGroups);
%  % stdRawODCountGroups  = reshape(avgDataset.IndivTwzrStdRawODCounts,  n, numGroups);
%  % 
%  % MeanODCountGroups = reshape(avgDataset.IndivTwzrMeanODCounts, n, numGroups);
%  % stdODCountGroups  = reshape(avgDataset.IndivTwzrStdODCounts,  n, numGroups);
%  % 
%  % % Mean of each group
%  % AllTwzrMeanRawODCounts = mean(MeanRawODCountGroups, 1)';
%  % AllTwzrMeanODCounts = mean(MeanODCountGroups, 1)';
%  % 
%  % % Propagated standard error of each mean
%  % AllTwzrStdRawODCounts = sqrt(sum(stdRawODCountGroups .^2, 1))' / n;
%  % AllTwzrStdODCounts = sqrt(sum(stdODCountGroups.^2, 1))' / n;
%  % 
%  % avgDataset.AllTwzrMeanRawODCounts   = AllTwzrMeanRawODCounts;
%  % avgDataset.AllTwzrStdRawODCounts    = AllTwzrStdRawODCounts;
%  % avgDataset.AllTwzrMeanODCounts      = AllTwzrMeanODCounts;
%  % avgDataset.AllTwzrStdODCounts       = AllTwzrStdODCounts;
%  
%  %% Average the mean counts across all tweezers with the same paramVal
%  %
%  % For each tweezer, the stored STD is the scatter of the individual
%  % measurements contributing to that tweezer's average.
%  %
%  % First convert each tweezer's STD to the standard error of its mean:
%  %
%  %       SEM_i = STD_i / sqrt(N_i)
%  %
%  % Then combine the tweezers using an unweighted mean.
%  %
%  % The uncertainty of the final mean contains two contributions:
%  %
%  %   1. Propagated uncertainty from the individual tweezer means:
%  %
%  %          propagatedError = sqrt(sum(SEM_i^2)) / N
%  %
%  %   2. Scatter between the individual tweezer means:
%  %
%  %          scatterError = std(tweezerMeans) / sqrt(N)
%  %
%  % The two contributions are combined in quadrature:
%  %
%  %          totalError = sqrt(propagatedError^2 + scatterError^2)
%  
%  
%  %% Number of parameter groups
%  
%  numGroups = avgDatasetCounter / numTweezers;
%  
%  
%  %% Reshape individual tweezer results into:
%  %
%  %       rows    = tweezers
%  %       columns = parameter groups
%  %
%  % This assumes the entries were created in the order:
%  %
%  %       parameter 1: tweezer 1, 2, 3, ...
%  %       parameter 2: tweezer 1, 2, 3, ...
%  %       etc.
%  
%  MeanRawODCountGroups = reshape( ...
%      avgDataset.IndivTwzrMeanRawODCounts, ...
%      numTweezers, ...
%      numGroups);
%  
%  StdRawODCountGroups = reshape( ...
%      avgDataset.IndivTwzrStdRawODCounts, ...
%      numTweezers, ...
%      numGroups);
%  
%  MeanODCountGroups = reshape( ...
%      avgDataset.IndivTwzrMeanODCounts, ...
%      numTweezers, ...
%      numGroups);
%  
%  StdODCountGroups = reshape( ...
%      avgDataset.IndivTwzrStdODCounts, ...
%      numTweezers, ...
%      numGroups);
%  
%  
%  %% Number of measurements contributing to each tweezer average
%  
%  NumAveragedGroups = reshape( ...
%      avgDataset.IndivTwzrNumAveraged, ...
%      numTweezers, ...
%      numGroups);
%  
%  
%  %% Convert individual-tweezer STD to SEM
%  
%  SemRawODCountGroups = ...
%      StdRawODCountGroups ./ sqrt(NumAveragedGroups);
%  
%  SemODCountGroups = ...
%      StdODCountGroups ./ sqrt(NumAveragedGroups);
%  
%  
%  %% Mean across tweezers
%  
%  AllTwzrMeanRawODCounts = ...
%      mean(MeanRawODCountGroups,1,'omitnan')';
%  
%  AllTwzrMeanODCounts = ...
%      mean(MeanODCountGroups,1,'omitnan')';
%  
%  
%  %% Calculate propagated uncertainty of the mean
%  %
%  % For an unweighted mean of N independent measurements:
%  %
%  %       propagatedError = sqrt(sum(SEM_i^2)) / N
%  %
%  % Ignore NaN values when determining which tweezers contribute.
%  
%  AllTwzrPropagatedRawError = nan(numGroups,1);
%  AllTwzrPropagatedODError  = nan(numGroups,1);
%  
%  for groupNum = 1:numGroups
%  
%      %% Raw counts
%  
%      semValues = SemRawODCountGroups(:,groupNum);
%  
%      valid = isfinite(semValues);
%  
%      numValid = sum(valid);
%  
%      if numValid > 0
%  
%          AllTwzrPropagatedRawError(groupNum) = ...
%              sqrt(sum(semValues(valid).^2)) / numValid;
%  
%      end
%  
%  
%      %% OD counts
%  
%      semValues = SemODCountGroups(:,groupNum);
%  
%      valid = isfinite(semValues);
%  
%      numValid = sum(valid);
%  
%      if numValid > 0
%  
%          AllTwzrPropagatedODError(groupNum) = ...
%              sqrt(sum(semValues(valid).^2)) / numValid;
%  
%      end
%  
%  end
%  
%  
%  %% Calculate scatter between tweezers
%  %
%  % The standard error associated with the observed spread of the
%  % individual tweezer means is:
%  %
%  %       scatterError = std(tweezerMeans) / sqrt(N)
%  
%  AllTwzrScatterRawError = nan(numGroups,1);
%  AllTwzrScatterODError  = nan(numGroups,1);
%  
%  for groupNum = 1:numGroups
%  
%      %% Raw counts
%  
%      values = MeanRawODCountGroups(:,groupNum);
%  
%      valid = isfinite(values);
%  
%      values = values(valid);
%  
%      numValid = numel(values);
%  
%      if numValid > 1
%  
%          AllTwzrScatterRawError(groupNum) = ...
%              std(values,0) / sqrt(numValid);
%  
%      elseif numValid == 1
%  
%          AllTwzrScatterRawError(groupNum) = 0;
%  
%      end
%  
%  
%      %% OD counts
%  
%      values = MeanODCountGroups(:,groupNum);
%  
%      valid = isfinite(values);
%  
%      values = values(valid);
%  
%      numValid = numel(values);
%  
%      if numValid > 1
%  
%          AllTwzrScatterODError(groupNum) = ...
%              std(values,0) / sqrt(numValid);
%  
%      elseif numValid == 1
%  
%          AllTwzrScatterODError(groupNum) = 0;
%  
%      end
%  
%  end
%  
%  
%  %% Combine propagated uncertainty and tweezer-to-tweezer scatter
%  %
%  %       totalError =
%  %           sqrt(propagatedError^2 + scatterError^2)
%  
%  AllTwzrStdRawODCounts = sqrt( ...
%      AllTwzrPropagatedRawError.^2 + ...
%      AllTwzrScatterRawError.^2);
%  
%  AllTwzrStdODCounts = sqrt( ...
%      AllTwzrPropagatedODError.^2 + ...
%      AllTwzrScatterODError.^2);
%  
%  
%  %% Store results
%  
%  avgDataset.AllTwzrMeanRawODCounts = ...
%      AllTwzrMeanRawODCounts;
%  
%  avgDataset.AllTwzrStdRawODCounts = ...
%      AllTwzrStdRawODCounts;
%  
%  avgDataset.AllTwzrMeanODCounts = ...
%      AllTwzrMeanODCounts;
%  
%  avgDataset.AllTwzrStdODCounts = ...
%      AllTwzrStdODCounts;
%  
%  
%  %% Also store the individual error contributions
%  %
%  % These aren't strictly necessary for plotting, but are extremely
%  % useful for debugging/understanding where the final error comes from.
%  
%  avgDataset.AllTwzrPropagatedRawError = ...
%      AllTwzrPropagatedRawError;
%  
%  avgDataset.AllTwzrScatterRawError = ...
%      AllTwzrScatterRawError;
%  
%  avgDataset.AllTwzrPropagatedODError = ...
%      AllTwzrPropagatedODError;
%  
%  avgDataset.AllTwzrScatterODError = ...
%      AllTwzrScatterODError;


fprintf(['Averaged OD images completed. ' ...
         'Created %d averaged entries.\n\n'], ...
    avgDatasetCounter);

fprintf('Averaged OD images completed. Created %d averaged entries.\n\n', avgDatasetCounter);

end