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

load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
numTweezers = size(tweezerROI.centersXY,1);

% From each individual file in indivDatset, Collect all values into a single list
% to scan over later when creating the properly averaged lists:
%   scanned parameter values (Or Dummy Scan file ID)
%   OD images
%   Total OD Counts
allParamVals    = [];
allTweezerNums  = [];
allImages       = {};
allSources      = {};
allODFiles      = {};
allTotODCounts  = [];
allTotODCountsImg1Raw = [];

allBasenameNums = [];
allImageNums    = [];

for basenameNum = 1:analyVar.numBasenamesAtom

    paramVals = get_average_scan_param_values(analyVar, indivDataset, basenameNum);

    for k = 1:indivDataset{basenameNum}.CounterAtom

        for tweezerNum = 1:numTweezers

            odFile = [analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                sprintf('_Tweezer%03d',tweezerNum) ...
                analyVar.ODimageFilename];

            if ~exist(odFile,'file')
                error('Missing OD image file:\n%s', odFile);
            end
    
            OD = dlmread(odFile);


            allParamVals(end+1,1)    = paramVals(k);
            allTweezerNums(end+1,1)  = tweezerNum;
            allImages{end+1,1}       = OD;
            allSources{end+1,1}      = char(indivDataset{basenameNum}.fileAtom(k));
            allODFiles{end+1,1}      = odFile;
            allBasenameNums(end+1,1) = basenameNum;
            allImageNums(end+1,1)    = k;
            allTotODCounts(end+1,1) = indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum);                   % PCA Bkg Subtracted Counts (May be normalized if toggled)
            allTotODCountsImg1Raw(end+1,1) = indivDataset{basenameNum}.OD_TotalCountsImg1Raw(k,tweezerNum);     % Raw Counts
        end
    end
end

%    % Group by parameter value
%    %Will group by common independent values (avg across multiple files) or
%    %group by all scans in a single file (avg a single dummy scan)
%    uniqueVals = unique(allParamVals);
%    
%    %% Initialize avgDataset
%    
%    % Any time there is a reference like: indivDataset{j}.whatever
%    % Instead use
%    %        src = avgDataset.sourceInfo{j};
%    %    
%    %        basenameNum = src.basenameNums(1);
%    %        k           = src.imageNums(1);
%    %        tweezerNum  = src.tweezerNums(1);
%    %        
%    %        value = indivDataset{basenameNum}.whatever(k);
%    
%    avgDataset = struct;
%    
%    avgDataset.avgDir       = avgDir;
%    avgDataset.numTweezers  = numTweezers;
%    avgDataset.TweezerROI   = tweezerROI;
%    
%    avgDataset.CounterAtom  = 0;
%    
%    avgDataset.paramVals    = [];
%    avgDataset.groupVals    = [];
%    avgDataset.tweezerNums  = [];
%    
%    avgDataset.avgODFiles   = {};
%    avgDataset.avgODImages  = {};
%    avgDataset.sourceFiles  = {};
%    avgDataset.sourceInfo   = {};
%    avgDataset.totODCounts  = {};
%    avgDataset.meanODCounts = [];
%    avgDataset.stdODCounts  = [];
%    avgDataset.numAveraged  = [];
%    
%    % Keep these for compatibility with older plotting/fitting functions
%    avgDataset.roiWin_Index = indivDataset{1}.roiWin_Index;
%    avgDataset.imagevcoAtom = [];
%    
%    %% Build averaged entries
%    avgCounter = 0;
%    
%    for p = 1:numel(uniqueVals)
%    
%        for tweezerNum = 1:numTweezers
%    
%            idx = find(allParamVals == uniqueVals(p) & ...
%                       allTweezerNums == tweezerNum);
%    
%            if isempty(idx)
%                continue;
%            end
%    
%            avgCounter = avgCounter + 1;
%    
%            imgStack = [];
%    
%            % for n = 1:numel(idx)
%            %     imgStack(:,:,n) = allImages{idx(n)};
%            % end
%    % 
%            % avgOD = mean(imgStack,3,'omitnan');
%            for n = 1:numel(idx)                  %%%%%%%%%%%%%%
%    
%                thisTweezer = allTweezerNums(idx(n));
%                %normFactor  = tweezerNormFactors(thisTweezer);
%            
%                imgStack(:,:,n) = allImages{idx(n)};% ./ normFactor;
%            
%            end
%            
%            avgOD = mean(imgStack,3,'omitnan');  %%%%%%%%%%%%
%    
%    
%            %%%%%%%%%%
%            rawCounts = allTotODCounts(idx);
%    
%            normCounts = rawCounts;
%            
%            for n = 1:numel(idx)
%                thisTweezer = allTweezerNums(idx(n));
%                normCounts(n) = rawCounts(n);% ./ tweezerNormFactors(thisTweezer);
%            end
%    
%            % Representative true parameter value before rounding/grouping
%            trueParamMean = mean(allParamVals(idx), 'omitnan');
%    
%            avgDataset.CounterAtom = avgCounter;
%    
%            avgDataset.paramVals(avgCounter,1)   = trueParamMean;
%            avgDataset.groupVals(avgCounter,1)   = uniqueVals(p);
%            avgDataset.tweezerNums(avgCounter,1) = tweezerNum;
%            avgDataset.imagevcoAtom(avgCounter,1) = trueParamMean;
%    
%            avgDataset.avgODImages{avgCounter,1} = avgOD;
%    
%            avgDataset.sourceFiles{avgCounter,1} = allSources(idx);
%            % avgDataset.totODCounts{avgCounter,1} = allTotODCounts(idx);
%    % 
%            % avgDataset.meanODCounts(avgCounter,1) = mean(allTotODCounts(idx), 'omitnan');
%            % avgDataset.stdODCounts(avgCounter,1)  = std(allTotODCounts(idx), 'omitnan');
%            avgDataset.numAveraged(avgCounter,1)  = numel(idx);
%    
%            avgDataset.rawODCounts{avgCounter,1}  = rawCounts;
%            avgDataset.normODCounts{avgCounter,1} = normCounts;
%            
%            avgDataset.totODCounts{avgCounter,1} = normCounts;
%            
%            avgDataset.meanRawODCounts(avgCounter,1) = mean(rawCounts, 'omitnan');
%            avgDataset.stdRawODCounts(avgCounter,1)  = std(rawCounts, 'omitnan');
%            
%            avgDataset.meanODCounts(avgCounter,1) = mean(normCounts, 'omitnan');
%            avgDataset.stdODCounts(avgCounter,1)  = std(normCounts, 'omitnan');
%    
%            % This is the important provenance block
%            sourceInfo = struct;
%            sourceInfo.basenameNums = allBasenameNums(idx);
%            sourceInfo.imageNums    = allImageNums(idx);
%            sourceInfo.tweezerNums  = allTweezerNums(idx);
%            sourceInfo.paramVals    = allParamVals(idx);
%            sourceInfo.groupVal     = uniqueVals(p);
%            sourceInfo.sourceFiles  = allSources(idx);
%            sourceInfo.odFiles      = allODFiles(idx);
%            %sourceInfo.totODCounts  = allTotODCounts(idx);
%            sourceInfo.rawODCounts = rawCounts; %%%%%%%%%
%            sourceInfo.normODCounts = normCounts;
%            sourceInfo.totODCounts = normCounts;
%            %sourceInfo.tweezerNormFactors = tweezerNormFactors(allTweezerNums(idx)); %%%%%%%%%%%
%    
%            avgDataset.sourceInfo{avgCounter,1} = sourceInfo;
%    
%            safeVal = regexprep(num2str(uniqueVals(p),'%.12g'), ...
%                '[^a-zA-Z0-9_\-\.]', '_');
%    
%            avgFile = fullfile(avgDataset.avgDir, ...
%                sprintf('AvgOD_%s_%s_Tweezer%03d.txt', ...
%                analyVar.avgScanParamField, safeVal, tweezerNum));
%    
%            dlmwrite(avgFile, avgOD, '\t');
%    
%            avgDataset.avgODFiles{avgCounter,1} = avgFile;
%    
%            fprintf('Avg entry %d: parameter %g, tweezer %d, using %d images\n', ...
%                avgCounter, trueParamMean, tweezerNum, numel(idx));
%    
%        end
%    end
%    
%    avgDataset.NormalizeTweezers = analyVar.NormalizeTweezers;
%    % avgDataset.tweezerNormFactors = tweezerNormFactors;
%    % avgDataset.tweezerMeanCounts = tweezerMeanCounts;

%% Group by scan parameter value

% Choose a tolerance much smaller than the spacing between scan points.
% Adjust this value if required for your parameter units.
paramTolerance = 1e-9;

if ~isvector(allParamVals) || ...
        ~isvector(allTweezerNums) || ...
        ~iscell(allImages)

    error(['allParamVals and allTweezerNums must be vectors, and ' ...
           'allImages must be a cell array.']);
end

allParamVals   = allParamVals(:);
allTweezerNums = allTweezerNums(:);
allTotODCounts = allTotODCounts(:);

numRecords = numel(allParamVals);

if numel(allTweezerNums) ~= numRecords || ...
        numel(allImages) ~= numRecords || ...
        numel(allTotODCounts) ~= numRecords

    error(['The flattened arrays do not contain the same number of ' ...
           'records:\n' ...
           '  allParamVals:    %d\n' ...
           '  allTweezerNums:  %d\n' ...
           '  allImages:       %d\n' ...
           '  allTotODCounts:  %d'], ...
        numel(allParamVals), ...
        numel(allTweezerNums), ...
        numel(allImages), ...
        numel(allTotODCounts));
end

%% Validate tweezer numbers

if any(~isfinite(allTweezerNums)) || ...
        any(allTweezerNums ~= round(allTweezerNums)) || ...
        any(allTweezerNums < 1) || ...
        any(allTweezerNums > numTweezers)

    error(['allTweezerNums contains an invalid tweezer number. ' ...
           'Expected integers from 1 through %d.'],numTweezers);
end

%% Create stable parameter-group values

groupedParamVals = ...
    round(allParamVals/paramTolerance)*paramTolerance;

uniqueVals = unique(groupedParamVals,'sorted');


%% Initialize avgDataset

avgDataset = struct;

avgDataset.avgDir       = avgDir;
avgDataset.numTweezers  = numTweezers;
avgDataset.TweezerROI   = tweezerROI;

avgDataset.CounterAtom  = 0;

avgDataset.paramVals    = [];
avgDataset.groupVals    = [];
avgDataset.tweezerNums  = [];

avgDataset.avgODFiles   = {};
avgDataset.avgODImages  = {};
avgDataset.sourceFiles  = {};
avgDataset.sourceInfo   = {};

avgDataset.rawODCounts       = {};
avgDataset.normODCounts      = {};
avgDataset.totODCounts       = {};
avgDataset.imageSumCounts    = {};

avgDataset.meanRawODCounts   = [];
avgDataset.stdRawODCounts    = [];

avgDataset.meanODCounts      = [];
avgDataset.stdODCounts       = [];

avgDataset.meanImageSumCounts = [];
avgDataset.stdImageSumCounts  = [];

avgDataset.numAveraged = [];

% Compatibility fields
avgDataset.roiWin_Index = indivDataset{1}.roiWin_Index;
avgDataset.imagevcoAtom = [];

% Record grouping information
avgDataset.paramTolerance   = paramTolerance;
avgDataset.uniqueGroupVals  = uniqueVals;


%% Build averaged parameter/tweezer entries

avgCounter = 0;

for p = 1:numel(uniqueVals)

    groupValue = uniqueVals(p);

    for tweezerNum = 1:numTweezers

        %% Find all records belonging to this exact group

        idx = find( ...
            groupedParamVals == groupValue & ...
            allTweezerNums == tweezerNum);

        if isempty(idx)

            warning(['No images found for parameter group %.12g, ' ...
                     'tweezer %d.'], ...
                groupValue,tweezerNum);

            continue;
        end

        %% Verify the selected records

        if any(allTweezerNums(idx) ~= tweezerNum)
            error(['Internal grouping error for parameter %.12g, ' ...
                   'tweezer %d.'], ...
                groupValue,tweezerNum);
        end

        selectedParamVals = allParamVals(idx);

        if any(abs(selectedParamVals-groupValue) > paramTolerance)
            error(['Parameter-grouping error for group %.12g, ' ...
                   'tweezer %d.'], ...
                groupValue,tweezerNum);
        end

        %% Verify image dimensions and create image stack

        firstImage = allImages{idx(1)};

        if isempty(firstImage) || ~ismatrix(firstImage)
            error(['Invalid image for source record %d, parameter %.12g, ' ...
                   'tweezer %d.'], ...
                idx(1),groupValue,tweezerNum);
        end

        imageRows = size(firstImage,1);
        imageCols = size(firstImage,2);

        imgStack = nan( ...
            imageRows, ...
            imageCols, ...
            numel(idx));

        imageSumCounts = nan(numel(idx),1);

        for n = 1:numel(idx)

            recordIndex = idx(n);
            thisImage = allImages{recordIndex};

            if ~isequal(size(thisImage),[imageRows imageCols])
                error(['Image-size mismatch for record %d. Expected ' ...
                       '%d-by-%d but found %d-by-%d.'], ...
                    recordIndex, ...
                    imageRows,imageCols, ...
                    size(thisImage,1),size(thisImage,2));
            end

            imgStack(:,:,n) = thisImage;

            % Counts obtained directly from the exact image being averaged
            imageSumCounts(n) = ...
                sum(thisImage(:),'omitnan');
        end

        %% Average images

        avgOD = mean(imgStack,3,'omitnan');

        %% Collect externally calculated counts

        rawCounts = allTotODCounts(idx);

        % Placeholder for optional normalization
        normCounts = rawCounts;

        %% Compare stored count values to sums of the actual images

        countDifference = ...
            rawCounts - imageSumCounts;

        finiteComparison = ...
            isfinite(rawCounts) & isfinite(imageSumCounts);

        if any(finiteComparison)

            comparisonScale = max( ...
                abs(rawCounts(finiteComparison)), ...
                abs(imageSumCounts(finiteComparison)));

            comparisonTolerance = ...
                1e-8*max(1,comparisonScale);

            badComparison = ...
                abs(countDifference(finiteComparison)) > ...
                comparisonTolerance;

            if any(badComparison)

                fprintf(['WARNING: image sums and allTotODCounts differ ' ...
                         'for parameter %.12g, tweezer %d.\n'], ...
                    groupValue,tweezerNum);

                fprintf(['  Mean allTotODCounts:       %.12g\n' ...
                         '  Mean direct image sum:     %.12g\n' ...
                         '  Mean difference:           %.12g\n'], ...
                    mean(rawCounts,'omitnan'), ...
                    mean(imageSumCounts,'omitnan'), ...
                    mean(countDifference,'omitnan'));
            end
        end

        %% Representative true parameter value

        trueParamMean = ...
            mean(selectedParamVals,'omitnan');

        %% Store averaged entry

        avgCounter = avgCounter + 1;

        avgDataset.CounterAtom = avgCounter;

        avgDataset.paramVals(avgCounter,1) = ...
            trueParamMean;

        avgDataset.groupVals(avgCounter,1) = ...
            groupValue;

        avgDataset.tweezerNums(avgCounter,1) = ...
            tweezerNum;

        % One entry for every averaged parameter/tweezer image
        avgDataset.imagevcoAtom(avgCounter,1) = ...
            trueParamMean;

        avgDataset.avgODImages{avgCounter,1} = ...
            avgOD;

        avgDataset.sourceFiles{avgCounter,1} = ...
            allSources(idx);

        avgDataset.numAveraged(avgCounter,1) = ...
            numel(idx);

        %% Store count information

        avgDataset.rawODCounts{avgCounter,1} = ...
            rawCounts;

        avgDataset.normODCounts{avgCounter,1} = ...
            normCounts;

        avgDataset.totODCounts{avgCounter,1} = ...
            normCounts;

        avgDataset.imageSumCounts{avgCounter,1} = ...
            imageSumCounts;

        avgDataset.meanRawODCounts(avgCounter,1) = ...
            mean(rawCounts,'omitnan');

        avgDataset.stdRawODCounts(avgCounter,1) = ...
            std(rawCounts,0,'omitnan');

        avgDataset.meanODCounts(avgCounter,1) = ...
            mean(normCounts,'omitnan');

        avgDataset.stdODCounts(avgCounter,1) = ...
            std(normCounts,0,'omitnan');

        avgDataset.meanImageSumCounts(avgCounter,1) = ...
            mean(imageSumCounts,'omitnan');

        avgDataset.stdImageSumCounts(avgCounter,1) = ...
            std(imageSumCounts,0,'omitnan');

        %% Store provenance

        sourceInfo = struct;

        sourceInfo.recordIndices = idx;

        sourceInfo.basenameNums = ...
            allBasenameNums(idx);

        sourceInfo.imageNums = ...
            allImageNums(idx);

        sourceInfo.tweezerNums = ...
            allTweezerNums(idx);

        sourceInfo.paramVals = ...
            allParamVals(idx);

        sourceInfo.groupedParamVals = ...
            groupedParamVals(idx);

        sourceInfo.groupVal = ...
            groupValue;

        sourceInfo.sourceFiles = ...
            allSources(idx);

        sourceInfo.odFiles = ...
            allODFiles(idx);

        sourceInfo.rawODCounts = ...
            rawCounts;

        sourceInfo.normODCounts = ...
            normCounts;

        sourceInfo.totODCounts = ...
            normCounts;

        sourceInfo.imageSumCounts = ...
            imageSumCounts;

        sourceInfo.countDifference = ...
            countDifference;

        avgDataset.sourceInfo{avgCounter,1} = ...
            sourceInfo;

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

        avgDataset.avgODFiles{avgCounter,1} = ...
            avgFile;

        fprintf(['Avg entry %d: parameter group %.12g, ' ...
                 'true mean %.12g, tweezer %d, using %d images\n'], ...
            avgCounter, ...
            groupValue, ...
            trueParamMean, ...
            tweezerNum, ...
            numel(idx));
    end
end

avgDataset.NormalizeTweezers = ...
    analyVar.NormalizeTweezers;

fprintf(['Averaged OD images completed. ' ...
         'Created %d averaged entries.\n\n'], ...
    avgCounter);

fprintf('Averaged OD images completed. Created %d averaged entries.\n\n', avgCounter);

end