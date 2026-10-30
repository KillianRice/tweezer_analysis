function imagefit_save_tweezer_spots(varargin)
%% Changed from imagefit_Backgrounds_PCA_V2 to just accomodate tweezer image analysis
%This program is designed to read several datafiles and corresponding background
%files, make ROI cuts based on tweezer spot location and save these images
%
% INPUTS:
%   varargin - variable input argument to allow passing of analysis
%              variables from Analysis_Runner. If not passed, the
%              program will call AnalysisData and get_indiv_batch_data itself.
%              It is important to follow the input construction below for
%              varargin to retrieve variable data from other programs.
%              -- first  argument - analyVar
%              -- second argument - indivDataset
%
% OUTPUTS:
%   none
%
% MISC:
%#ok<*PFBNS> - suppress all instances of 'this variables is indexed but not
%              sliced'. Passing large arrays may incur unnecessary communication
%              overhead to the workers but we are only reading from the
%              analyVar and indivDataset structures and not changing them.
%              If there is a performance issue look into the function
%              WorkerObjWrapper
%
% NOTES:
%   10.13.2026 - Changed code from imagefit_Background_PCA_V2 for
%   readability and processing spceficially for tweezers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Workspace specific options
warning('off','all') % Suppress parfor warnings in command window

%% Load variables and file data
if nargin == 0 % If run without arguments
    analyVar     = AnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
    avgBkgLibrary = create_AvgBkg_Library(analyVar, indivDataset);
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
    avgBkgLibrary = varargin{3}; % avg bkg image library is third
end

% %% Initialize parallel workers
% if isempty(gcp('nocreate'))
%     parpool('local', 2);
% end

 %% Load the parameters of the ROIs made by
 %%% check_pnts_tweezer.m
 roiData = load(fullfile(analyVar.analyOutDir, 'tweezerROI.mat'), 'tweezerROI');
 tweezerROI = roiData.tweezerROI;
 numTweezers = size(tweezerROI.centersXY,1);

%% Loop through each batch file listed in analyVar.basenamevectorAtom
for basenameNum = 1:analyVar.numBasenamesAtom
    % this will keep track of all the files analyzed in all the batches
    fprintf('\nSaving Image batch file %g of %g\n',basenameNum,analyVar.numBasenamesAtom)
    
    
%% Loop through each image and find background state
%%%%%%%-----------------------------------%%%%%%%%%%
    % this loop cuts ROIs for all atom and bkg files in this dataset
    for k = 1:indivDataset{basenameNum}.CounterAtom  

%% Create ROIs of raw and subtracted (OD) image
%%%%%%%-----------------------------------%%%%%%%%%%
        % Retrieve cloud matrix for single image
        roiImageAtom = reshape(indivDataset{basenameNum}.AtomsCloud(:,k),[1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1));
        if analyVar.avgBkgImg
            % idx = find([avgBkgLibrary.scanID] == analyVar.meanListVar(basenameNum) & ...
            %     [avgBkgLibrary.scanParameter] == indivDataset{basenameNum}.imagevcoAtom(k), 1);
            idx = find([avgBkgLibrary.scanID] == analyVar.meanListVar(basenameNum));
            if ~isempty(idx)
                avgBkgImage = avgBkgLibrary(idx).avgImageBkg;
            else
                disp('No matching averaged image found.');
            end
            roiImageBack = reshape(avgBkgImage(:),[1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1));
        else
            roiImageBack = reshape(indivDataset{basenameNum}.BackgroundCloud(:,k),[1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1));
        end

        % Bin and trim atoms images
        cutImageCell = TrimAndBin(analyVar,mat2cell([roiImageAtom;roiImageBack],[1 1]*size(roiImageAtom,1),size(roiImageAtom,2)));
        [cutImageAtom, cutImageBack] = cutImageCell{:};
        
        %% Subtract a fitted background
        if analyVar.UseImages_Fluorescence == 1
            disp('Fluorescence Images: subtracting background')
        
            OD_Image_Single = cutImageAtom - cutImageBack; 
        else
            err('Turn Fluorescence ON');

        end    

        %% When Using Tweezers: Multiple ROIs in Single Image
        % Note each tweezer ROI background is from the avg bkg image (not
        % raw)
        if analyVar.UseTweezer == 1

            %% Create 3D matrix for saving all ROIs
            roiSize = 2*tweezerROI.roiHalfWidthPix + 1;
            
            RawAtomROIs = zeros(roiSize, roiSize, numTweezers, 'like', cutImageAtom);
            RawBkgROIs  = zeros(roiSize, roiSize, numTweezers, 'like', cutImageBack);
            ODROIs      = zeros(roiSize, roiSize, numTweezers, 'like', OD_Image_Single);
        
            %%% For each image within a file, saving multiple images from
            %%% the ROI cuts
            for tweezerNum = 1:numTweezers

                x0 = tweezerROI.centersXY(tweezerNum,1);  %row
                y0 = tweezerROI.centersXY(tweezerNum,2); %column
                r  = tweezerROI.roiHalfWidthPix;

                rows = x0-r:x0+r;
                cols = y0-r:y0+r;
        
                %% Get ROI cut from image - bkg image:  Save individual txt files (uncomment if want txt)
                % OD_Tweezer = OD_Image_Single( rows, cols);
                % RawAtomImage_Tweezer = cutImageAtom( rows, cols);
                % RawBkgImage_Tweezer = cutImageBack( rows, cols);
                % 
                % if isempty(OD_Tweezer)
                %     error('Empty tweezer OD crop for tweezer %d.', tweezerNum);
                % end
                % 
                % % Save tweezer ROI bkg subtracted image
                % saveName = [analyVar.analyOutDir ...
                %     char(indivDataset{basenameNum}.fileAtom(k)) ...
                %     sprintf('_Tweezer%03d',tweezerNum) ...
                %     analyVar.ODimageFilename];
                %     
                % writematrix(OD_Tweezer, saveName, 'Delimiter', '\t');
                % % Save tweezer ROI raw atom image (image 1)
                % saveName = [analyVar.analyOutDir ...
                %     char(indivDataset{basenameNum}.fileAtom(k)) ...
                %     sprintf('_TweezerImg1Raw%03d',tweezerNum) ...
                %     analyVar.ODimageFilename];
                % 
                % writematrix(RawAtomImage_Tweezer, saveName, 'Delimiter', '\t');
                % % Save tweezer ROI raw bkg image (image 2)
                % saveName = [analyVar.analyOutDir ...
                %     char(indivDataset{basenameNum}.fileBack(k)) ...
                %     sprintf('_TweezerImg2Raw%03d',tweezerNum) ...
                %     analyVar.ODimageFilename];
                % 
                % writematrix(RawBkgImage_Tweezer, saveName, 'Delimiter', '\t');

                % Crop images, save into 3D matrix
                ODROIs(:,:,tweezerNum) = OD_Image_Single(rows, cols);
                RawAtomROIs(:,:,tweezerNum) = cutImageAtom(rows, cols);
                RawBkgROIs(:,:,tweezerNum) = cutImageBack(rows, cols);

            end
            % Save all tweezers for this scan in ONE file
            roiFile = fullfile(analyVar.analyOutDir, ...
                [char(indivDataset{basenameNum}.fileAtom(k)) '_ROIs.mat']);
            
            % mat file holds 3 3D arrays each indexed by tweezer
            save(roiFile, 'ODROIs', 'RawAtomROIs', 'RawBkgROIs', '-v7');
        else

        
            %% Save OD Single Image 
            dlmwrite([analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                analyVar.ODimageFilename], OD_Image_Single, '\t'); % (LEGACY)
        end
        

        
    end %loop through files in dataset (one batchfile)
end %end loop through batches

%% Cleanup Parallel workspace
warning('on','all') %% Reenable warnings
% delete(gcp('nocreate'));

%% Wrap Up
fclose('all'); % Close any file handles which may be open
fprintf('Saving is completed.\n\n')
end