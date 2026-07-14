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

% From each individual file in indivDatset, Collect all:
% scanned parameter values (Or Dummy Scan file ID)
% OD images
% Total OD Counts
allParamVals    = [];
allTweezerNums  = [];
allImages       = {};
allSources      = {};
allODFiles      = {};
allTotODCounts  = [];

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
            allTotODCounts(end+1,1) = indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum);
        end
    end
end

%% Tweezer normalization factors
if isfield(analyVar,'NormalizeTweezers') && analyVar.NormalizeTweezers == 1

    tweezerMeanCounts = zeros(numTweezers,1);

    for tweezerNum = 1:numTweezers
        idxT = allTweezerNums == tweezerNum;
        tweezerMeanCounts(tweezerNum) = mean(allTotODCounts(idxT), 'omitnan');
    end

    globalMeanCounts = mean(tweezerMeanCounts, 'omitnan');

    tweezerNormFactors = tweezerMeanCounts ./ globalMeanCounts;

    % Avoid divide-by-zero
    tweezerNormFactors(tweezerNormFactors == 0 | isnan(tweezerNormFactors)) = 1;

else

    tweezerNormFactors = ones(numTweezers,1);

end

% Group by parameter value
%Will group by common independent values (avg across multiple files) or
%group by all scans in a single file (avg a single dummy scan)
uniqueVals = unique(allParamVals);

%% Initialize avgDataset

% Any time there is a reference like: indivDataset{j}.whatever
% Instead use
%        src = avgDataset.sourceInfo{j};
%    
%        basenameNum = src.basenameNums(1);
%        k           = src.imageNums(1);
%        tweezerNum  = src.tweezerNums(1);
%        
%        value = indivDataset{basenameNum}.whatever(k);

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
avgDataset.totODCounts  = {};
avgDataset.meanODCounts = [];
avgDataset.stdODCounts  = [];
avgDataset.numAveraged  = [];

% Keep these for compatibility with older plotting/fitting functions
avgDataset.roiWin_Index = indivDataset{1}.roiWin_Index;
avgDataset.imagevcoAtom = [];

%% Build averaged entries
avgCounter = 0;

for p = 1:numel(uniqueVals)

    for tweezerNum = 1:numTweezers

        idx = find(allParamVals == uniqueVals(p) & ...
                   allTweezerNums == tweezerNum);

        if isempty(idx)
            continue;
        end

        avgCounter = avgCounter + 1;

        imgStack = [];

        % for n = 1:numel(idx)
        %     imgStack(:,:,n) = allImages{idx(n)};
        % end
% 
        % avgOD = mean(imgStack,3,'omitnan');
        for n = 1:numel(idx)                  %%%%%%%%%%%%%%

            thisTweezer = allTweezerNums(idx(n));
            normFactor  = tweezerNormFactors(thisTweezer);
        
            imgStack(:,:,n) = allImages{idx(n)} ./ normFactor;
        
        end
        
        avgOD = mean(imgStack,3,'omitnan');  %%%%%%%%%%%%


        %%%%%%%%%%
        rawCounts = allTotODCounts(idx);

        normCounts = rawCounts;
        
        for n = 1:numel(idx)
            thisTweezer = allTweezerNums(idx(n));
            normCounts(n) = rawCounts(n) ./ tweezerNormFactors(thisTweezer);
        end

        % Representative true parameter value before rounding/grouping
        trueParamMean = mean(allParamVals(idx), 'omitnan');

        avgDataset.CounterAtom = avgCounter;

        avgDataset.paramVals(avgCounter,1)   = trueParamMean;
        avgDataset.groupVals(avgCounter,1)   = uniqueVals(p);
        avgDataset.tweezerNums(avgCounter,1) = tweezerNum;
        avgDataset.imagevcoAtom(avgCounter,1) = trueParamMean;

        avgDataset.avgODImages{avgCounter,1} = avgOD;

        avgDataset.sourceFiles{avgCounter,1} = allSources(idx);
        % avgDataset.totODCounts{avgCounter,1} = allTotODCounts(idx);
% 
        % avgDataset.meanODCounts(avgCounter,1) = mean(allTotODCounts(idx), 'omitnan');
        % avgDataset.stdODCounts(avgCounter,1)  = std(allTotODCounts(idx), 'omitnan');
        avgDataset.numAveraged(avgCounter,1)  = numel(idx);

        avgDataset.rawODCounts{avgCounter,1}  = rawCounts;
        avgDataset.normODCounts{avgCounter,1} = normCounts;
        
        avgDataset.totODCounts{avgCounter,1} = normCounts;
        
        avgDataset.meanRawODCounts(avgCounter,1) = mean(rawCounts, 'omitnan');
        avgDataset.stdRawODCounts(avgCounter,1)  = std(rawCounts, 'omitnan');
        
        avgDataset.meanODCounts(avgCounter,1) = mean(normCounts, 'omitnan');
        avgDataset.stdODCounts(avgCounter,1)  = std(normCounts, 'omitnan');

        % This is the important provenance block
        sourceInfo = struct;
        sourceInfo.basenameNums = allBasenameNums(idx);
        sourceInfo.imageNums    = allImageNums(idx);
        sourceInfo.tweezerNums  = allTweezerNums(idx);
        sourceInfo.paramVals    = allParamVals(idx);
        sourceInfo.groupVal     = uniqueVals(p);
        sourceInfo.sourceFiles  = allSources(idx);
        sourceInfo.odFiles      = allODFiles(idx);
        %sourceInfo.totODCounts  = allTotODCounts(idx);
        sourceInfo.rawODCounts = rawCounts; %%%%%%%%%
        sourceInfo.normODCounts = normCounts;
        sourceInfo.totODCounts = normCounts;
        sourceInfo.tweezerNormFactors = tweezerNormFactors(allTweezerNums(idx)); %%%%%%%%%%%

        avgDataset.sourceInfo{avgCounter,1} = sourceInfo;

        safeVal = regexprep(num2str(uniqueVals(p),'%.12g'), ...
            '[^a-zA-Z0-9_\-\.]', '_');

        avgFile = fullfile(avgDataset.avgDir, ...
            sprintf('AvgOD_%s_%s_Tweezer%03d.txt', ...
            analyVar.avgScanParamField, safeVal, tweezerNum));

        dlmwrite(avgFile, avgOD, '\t');

        avgDataset.avgODFiles{avgCounter,1} = avgFile;

        fprintf('Avg entry %d: parameter %g, tweezer %d, using %d images\n', ...
            avgCounter, trueParamMean, tweezerNum, numel(idx));

    end
end

avgDataset.NormalizeTweezers = analyVar.NormalizeTweezers;
avgDataset.tweezerNormFactors = tweezerNormFactors;
avgDataset.tweezerMeanCounts = tweezerMeanCounts;

fprintf('Averaged OD images completed. Created %d averaged entries.\n\n', avgCounter);

end

% avgDataset = struct;
% avgDataset.avgDir = avgDir;
% avgDataset.roiWin_Index = indivDataset{1}.roiWin_Index;
% avgDataset.imagevcoAtom = uniqueVals;
% avgDataset.CounterAtom = numel(uniqueVals);
% avgDataset.avgODFiles = cell(numel(uniqueVals),1);
% avgDataset.sourceFiles = cell(numel(uniqueVals),1);
% avgDataset.totODCounts = cell(numel(uniqueVals),1);
% 
% avgCounter = 0;
% 
% for p = 1:numel(uniqueVals)
% 
%     for tweezerNum = 1:numTweezers
% 
%         idx = find(allParamVals == uniqueVals(p) & ...
%                    allTweezerNums == tweezerNum);
% 
%         if isempty(idx)
%             continue;
%         end
% 
%         avgCounter = avgCounter + 1;
% 
%         imgStack = [];
% 
%         for n = 1:numel(idx)
%             imgStack(:,:,n) = allImages{idx(n)};
%         end
% 
%         avgOD = mean(imgStack,3,'omitnan');
% 
%         avgDataset.paramVals(avgCounter,1) = uniqueVals(p);
%         avgDataset.tweezerNums(avgCounter,1) = tweezerNum;
%         avgDataset.avgODImages{avgCounter} = avgOD;
% 
%         safeVal = regexprep(num2str(uniqueVals(p),'%.12g'), ...
%             '[^a-zA-Z0-9_\-\.]', '_');
% 
%         avgFile = fullfile(avgDataset.avgDir, ...
%             sprintf('AvgOD_%s_%s_Tweezer%03d.txt', ...
%             analyVar.avgScanParamField, safeVal, tweezerNum));
% 
%         dlmwrite(avgFile, avgOD, '\t');
% 
%         avgDataset.avgODFiles{avgCounter} = avgFile;
%     end
% 
%     fprintf('Averaged parameter value %g using %d images\n', uniqueVals(p), numel(idx));
% end
% 
% 
% %save(fullfile(avgDir,'avgDataset.mat'),'avgDataset');
% 
% avgDataset.CounterAtom = avgCounter;
% avgDataset.numTweezers = numTweezers;
% avgDataset.TweezerROI = tweezerROI;
% 
% fprintf('Averaged OD images completed.\n\n');
% 
% end