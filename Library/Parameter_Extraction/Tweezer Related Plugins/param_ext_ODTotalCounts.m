function indivDataset = param_ext_ODTotalCounts(analyVar,indivDataset)
%
% Each file now has
%
% indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum)
% indivDataset{basenameNum}.TotalCountsImg1Raw(k,tweezerNum)
% indivDataset{basenameNum}.TotalCountsRawBkg(k,tweezerNum)   
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

fprintf('\nCounting Tweezer ROIs...\n');

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
                %% Read saved OD image for each tweezer & Bkg
                atomFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileAtom(k)) ...
                    sprintf('_Tweezer%03d',tweezerNum) ...
                    analyVar.ODimageFilename];
    
                if ~exist(atomFile,'file')
                    error('Missing tweezer OD image:\n%s', atomFile);
                end
    
                %OD_Image = dlmread(odFile);
                atom_Image = readmatrix(atomFile, 'FileType', 'text');

                if isempty(atom_Image)
                    error('OD counter found an empty OD image file:\n%s', atomFile);
                end

                %% Read saved Bkg image for each tweezer (Image2)
                bkgFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileBack(k)) ...
                    sprintf('_TweezerImg2Raw%03d',tweezerNum) ...
                    analyVar.ODimageFilename];
    
                if ~exist(bkgFile,'file')
                    error('Missing tweezer OD image:\n%s', bkgFile);
                end
    
                bkgImage = readmatrix(bkgFile, 'FileType', 'text');

                if isempty(bkgImage)
                    error('OD counter found an empty OD image file:\n%s', bkgFile);
                end

                %% Read saved Raw ROI Cut Atom Image (Image 1)
                atomFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileAtom(k)) ...
                    sprintf('_TweezerImg1Raw%03d',tweezerNum) ...
                    analyVar.ODimageFilename];
    
                if ~exist(atomFile,'file')
                    error('Missing tweezer OD image:\n%s', atomFile);
                end
    
                atom_Image = readmatrix(atomFile, 'FileType', 'text');

                if isempty(atom_Image)
                    error('OD counter found an empty OD image file:\n%s', atomFile);
                end

                % Save the full image to indivDataset
                indivDataset{basenameNum}.backSubtractedTwzImg{k,tweezerNum} = ...
                    atom_Image - bkgImage;

                % Save sum of the x-y cuts to the indivDataset
                indivDataset{basenameNum}.backSubtractedTwzImgXSum{k,tweezerNum} = ...
                    sum(atom_Image - bkgImage, 1).';
                indivDataset{basenameNum}.backSubtractedTwzImgYSum{k,tweezerNum} = ...
                    sum(atom_Image - bkgImage, 2);
                
                %% Integrated counts: Save atom-bkg = OD and also the raw bkg and atom counts
                %OD
                indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum) = ...
                    sum(atom_Image(:),'omitnan') - sum(bkgImage(:),'omitnan');
                %Bkg
                indivDataset{basenameNum}.TotalCountsRawBkg(k,tweezerNum) = ...
                    sum(bkgImage(:),'omitnan');
                
                %Atom
                indivDataset{basenameNum}.TotalCountsImg1Raw(k,tweezerNum) = ...
                    sum(atom_Image(:),'omitnan');

                %% Check if there was atom or no atom (using a given threshold)
                if indivDataset{basenameNum}.OD_TotalCounts(k,tweezerNum) >= analyVar.atomCountSignal
                    indivDataset{basenameNum}.isThereAtom(k,tweezerNum) = 1;
                else
                    indivDataset{basenameNum}.isThereAtom(k,tweezerNum) = 0;
                end
             end
        else
            %% Read saved OD image
            atomFile = [analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                analyVar.ODimageFilename];
    
            if ~exist(atomFile,'file')
                error('Cannot find OD image:\n%s',atomFile);
            end
    
            atom_Image = dlmread(atomFile);
    
            %% Integrated counts
            indivDataset{basenameNum}.OD_TotalCounts(k) = ...
                sum(atom_Image(:),'omitnan');

        end
    end

    for tweezerNum = 1:numTweezers
        indivDataset{basenameNum}.isThereAtomPercentage(tweezerNum) = ... 
                mean(indivDataset{basenameNum}.isThereAtom(:,tweezerNum));
    end
end

% Sum the tweezer image arrays and average
if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1

    fprintf('\nBuilding 1D Image Summed Arrays...\n');

    for basenameNum = 1:analyVar.numBasenamesAtom
        ImgXSum = zeros (size(indivDataset{basenameNum}.backSubtractedTwzImgYSum{1,tweezerNum},1),1); %% Create Pixel length Column Array
        ImgYSum = zeros (size(indivDataset{basenameNum}.backSubtractedTwzImgXSum{1,tweezerNum},1),1); %% Create Pixel length Column Array
        
         for k = 1:indivDataset{basenameNum}.CounterAtom
             for tweezerNum = 1:numTweezers
                 ImgXSum =  ImgXSum + indivDataset{basenameNum}.backSubtractedTwzImgYSum{k,tweezerNum};
                 ImgYSum =  ImgYSum + indivDataset{basenameNum}.backSubtractedTwzImgXSum{k,tweezerNum};
             end
         end
            %indivDataset{basenameNum}.backSubtractedTwzImgYSum{k,tweezerNum},1
         indivDataset{basenameNum}.ImgXSum = ImgXSum./ indivDataset{basenameNum}.CounterAtom;
         indivDataset{basenameNum}.ImgYSum = ImgYSum./ indivDataset{basenameNum}.CounterAtom;
         indivDataset{basenameNum}.ImgPixelLengthList = 1:size(ImgXSum,1);
    end
    
end 


%% Normalize the ODTotalCounts
if isfield(analyVar,'NormalizeTweezers') && analyVar.NormalizeTweezers == 1

    fprintf('\nBuilding Normalized Counts...\n');
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