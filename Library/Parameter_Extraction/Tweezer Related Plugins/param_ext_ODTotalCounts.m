function indivDataset = param_ext_ODTotalCounts(analyVar,indivDataset)
%
% Each file now has
%
% indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum)
% indivDataset{basenameNum}.OD_TotalCountsImg1Raw(k,tweezerNum)
% or
% indivDataset{basenameNum}.OD_TotalCounts(k)
%
% This will grab the Total OD counts within an image and store that within
% indivDataset to be grabbed later by imagefit_BuildAveragedScans when
% creating the correct collection of tweezer spots and scanned parameters

%% When Using Tweezers: Multiple ROIs in Single Image
if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
    %%% Grab Tweezer info
    load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
    numTweezers = size(tweezerROI.centersXY,1);
end

%% Loop over scans
for basenameNum = 1:analyVar.numBasenamesAtom

    % Preallocate
    %%% If using tweezer: preallocate Number of Images * Number of Tweezer
    if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
        indivDataset{basenameNum}.OD_TotalCounts = ...
        zeros(indivDataset{basenameNum}.CounterAtom, numTweezers);
    else
        indivDataset{basenameNum}.OD_TotalCounts = ...
            zeros(1, indivDataset{basenameNum}.CounterAtom);
    end

    for k = 1:indivDataset{basenameNum}.CounterAtom
        %%% If using tweezer:
        if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
             for tweezerNum = 1:numTweezers
                %% Read saved OD image for each tweezer
                odFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileAtom(k)) ...
                    sprintf('_Tweezer%03d',tweezerNum) ...
                    analyVar.ODimageFilename];
    
                if ~exist(odFile,'file')
                    error('Missing tweezer OD image:\n%s', odFile);
                end
    
                %OD_Image = dlmread(odFile);
                OD_Image = readmatrix(odFile, 'FileType', 'text');

                if isempty(OD_Image)
                    error('OD counter found an empty OD image file:\n%s', odFile);
                end
                
                 %% Integrated counts
                indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum) = ...
                    sum(OD_Image(:),'omitnan');

                %% Read saved Raw ROI Cut Atom Image (Image 1)
                odFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileAtom(k)) ...
                    sprintf('_TweezerImg1Raw%03d',tweezerNum) ...
                    analyVar.ODimageFilename];
    
                if ~exist(odFile,'file')
                    error('Missing tweezer OD image:\n%s', odFile);
                end
    
                %OD_Image = dlmread(odFile);
                OD_Image = readmatrix(odFile, 'FileType', 'text');

                if isempty(OD_Image)
                    error('OD counter found an empty OD image file:\n%s', odFile);
                end
                
                 %% Integrated counts
                indivDataset{basenameNum}.OD_TotalCountsImg1Raw(k,tweezerNum) = ...
                    sum(OD_Image(:),'omitnan');
             end
        else
            %% Read saved OD image
            odFile = [analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                analyVar.ODimageFilename];
    
            if ~exist(odFile,'file')
                error('Cannot find OD image:\n%s',odFile);
            end
    
            OD_Image = dlmread(odFile);
    
            %% Integrated counts
            indivDataset{basenameNum}.OD_TotalCounts(k) = ...
                sum(OD_Image(:),'omitnan');

        end
    end
end


%% Normalize the ODTotalCounts
if isfield(analyVar,'NormalizeTweezers') && analyVar.NormalizeTweezers == 1

    % first gather all the OD counts from the subtracted background ROI Tweezer
    % cuts.
    % Also keep track of which tweezer each entry comes from
    allTotODCounts  = [];
    allTweezerNums  = [];
    % loop through files
    for basenameNum = 1:analyVar.numBasenamesAtom
        %loop through images in a single file
        for k = 1:indivDataset{basenameNum}.CounterAtom
            % loop through tweezers within a single image
            for tweezerNum = 1:numTweezers
                allTotODCounts(end+1,1) = indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum);
                allTweezerNums(end+1,1)  = tweezerNum;
            end
        end
    end
    tweezerMeanCounts = zeros(numTweezers,1);

    % Find all indexes of a certain tweezer (such as top left is tweezer 1)
    for tweezerNum = 1:numTweezers
        idxT = allTweezerNums == tweezerNum;
        tweezerMeanCounts(tweezerNum) = mean(allTotODCounts(idxT), 'omitnan');
    end
    
    % Each Tweezer ROI will have its counts divided by the mean of all the
    % data collected for that ROI
    tweezerNormFactors = tweezerMeanCounts;

    % Avoid divide-by-zero
    tweezerNormFactors(tweezerNormFactors == 0 | isnan(tweezerNormFactors)) = 1;


    % Finally loop through the indivDataset again and normalize the counts
    for basenameNum = 1:analyVar.numBasenamesAtom
    
        for k = 1:indivDataset{basenameNum}.CounterAtom
                
            for tweezerNum = 1:numTweezers
                indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum) = indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum) ./ tweezerNormFactors(tweezerNum);
            end
        end
    end


end


end