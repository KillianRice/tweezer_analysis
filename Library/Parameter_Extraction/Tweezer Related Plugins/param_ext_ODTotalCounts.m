function indivDataset = param_ext_ODTotalCounts(analyVar,indivDataset)
%
% Each file now has
%
% indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum)
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

end