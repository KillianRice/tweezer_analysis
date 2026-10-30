function avgBkgLibrary = create_AvgBkg_Library(varargin)
%% BUILDAVGBKGLIBRARY
% Groups and averages background images by scan ID and scan parameter.
%
% INPUTS:
%   indivDataset : 1D cell array of structs, one per file
%   analyVar     : Analysis variables struct
%
% OUTPUTS:
%   avgBkgLibrary : Struct array containing averaged images and metadata
%
% EXPECTED INPUT FIELDS:
%   indivDataset{i}.rawBackImages
%   indivDataset{i}.CounterAtom
%   indivDataset{i}.imagevcoAtom
%   analyVar.meanListVar
%   analyVar.numBasenamesAtoms
%   analyVar.analyOutDir


%% Load variables and file data
if nargin == 0 % If run without arguments
    analyVar     = AnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
end


fprintf('\nCollecting Image Backgrounds...\n');

    %% 1. Determine number of files

    numFiles = analyVar.numBasenamesAtom;

    %% 2. Collect all images and their metadata
    
    allImagesAtom = {};
    allImagesBkg = {};
    allScanIDs = [];
    allScanParams = [];
    allSourceFiles = {};
    allSourceIndices = [];

    for i = 1:numFiles

        numImages = indivDataset{i}.CounterAtom;

        rawImagesAtom = indivDataset{i}.AtomsCloud;
        rawImagesBkg = indivDataset{i}.BackgroundCloud;
        scanParams = indivDataset{i}.imagevcoAtom;
        scanID = analyVar.meanListVar(i);
        % if dummy scan then all images in a file should be averaged to
        % each other. So just give the same scan param as the scan ID for
        % easy sorting in below machinery
        if analyVar.dummyScan
            for j = 1:length(scanParams)
                scanParams(j) = scanID;
            end
        end

        % Check that the file contains the expected number of images.
        if size(rawImagesBkg, 2) < numImages
            error('File %d has fewer image columns than CounterAtom.', i);
        end

        if numel(scanParams) < numImages
            error('File %d has fewer scan parameters than CounterAtom.', i);
        end

        % Determine a readable source filename.
        if isfield(indivDataset{i}, 'fileAtom') && ...
                ~isempty(indivDataset{i}.fileAtom)
            sourceFile = char(indivDataset{i}.fileAtom(1));
        elseif isfield(indivDataset{i}, 'fileBack') && ...
                ~isempty(indivDataset{i}.fileBack)
            sourceFile = char(indivDataset{i}.fileBack(1));
        else
            sourceFile = sprintf('File_%d', i);
        end

        % Store every individual image with its metadata.
        for k = 1:numImages
            
            allImagesAtom{end+1} = rawImagesAtom(:, k);
            allImagesBkg{end+1} = rawImagesBkg(:, k);

            allScanIDs(end+1, 1) = scanID;
            allScanParams(end+1, 1) = scanParams(k);

            allSourceFiles{end+1, 1} = sourceFile;
            allSourceIndices(end+1, 1) = k;

        end
    end

    if isempty(allImagesBkg)
        error('No background images were found.');
    end

    %% 3. Determine image dimensions

    numPixels = numel(allImagesBkg{1});

    for k = 1:numel(allImagesBkg)
        if numel(allImagesBkg{k}) ~= numPixels
            error('All background images must have the same dimensions.');
        end
    end

 

    sideLength = sqrt(numPixels);   
    if abs(sideLength - round(sideLength)) < 1e-10
        imageSize = [round(sideLength), round(sideLength)];
    else
        error(['Cannot infer image dimensions from %d pixels. ', ...
            'Pass imageSize, e.g. [200 200].'], numPixels);
    end

    %% 4. Group by scan ID and scan parameter

    % Each row is one unique (scan ID, scan parameter) combination.
    groupKeys = unique([allScanIDs, allScanParams], 'rows', 'stable');

    numGroups = size(groupKeys, 1);

    avgBkgLibrary = struct( ...
        'scanID', [], ...
        'scanParameter', [], ...
        'avgImageBkg', [], ...
        'avgImage2DBkg', [], ...
        'numImages', [], ...
        'sourceFiles', {{}}, ...
        'sourceImageIndices', [], ...
        'outputFiles', [], ...
        'avgImageSubtracted', []);

    for g = 1:numGroups

        thisScanID = groupKeys(g, 1);
        thisParam = groupKeys(g, 2);

        % Find all images matching BOTH values.
        matchIdx = ...
            (allScanIDs == thisScanID) & ...
            (allScanParams == thisParam);

        imageIdx = find(matchIdx);

        % Combine matching image vectors into a matrix.
        imageStackAtom = zeros(numPixels, numel(imageIdx));
        imageStackBkg = zeros(numPixels, numel(imageIdx));

        for j = 1:numel(imageIdx)
            imageStackAtom(:, j) = allImagesAtom{imageIdx(j)};
            imageStackBkg(:, j) = allImagesBkg{imageIdx(j)};
        end

        % Average across images (columns).
        meanImageAtom = mean(imageStackAtom, 2);
        meanImageBkg = mean(imageStackBkg, 2);
        stdAvgImageAtom = std(imageStackAtom, 0, 2);
        stdAvgImageBkg = std(imageStackBkg, 0, 2);

        % Populate the output struct.
        avgBkgLibrary(g).scanID = thisScanID;
        avgBkgLibrary(g).scanParameter = thisParam;
        
        avgBkgLibrary(g).avgImageAtom = meanImageAtom;
        avgBkgLibrary(g).stdAvgImageAtom = stdAvgImageAtom;
        avgBkgLibrary(g).avgImage2DAtom = reshape(meanImageAtom, imageSize);
        
        avgBkgLibrary(g).avgImageBkg = meanImageBkg;
        avgBkgLibrary(g).stdAvgImageBkg = stdAvgImageBkg;
        avgBkgLibrary(g).avgImage2DBkg = reshape(meanImageBkg, imageSize);
        
        avgBkgLibrary(g).numImages = numel(imageIdx);

        avgBkgLibrary(g).sourceFiles = ...
            unique(allSourceFiles(matchIdx), 'stable');

        avgBkgLibrary(g).sourceImageIndices = ...
            allSourceIndices(matchIdx);

        differenceImg = meanImageAtom - meanImageBkg;
        difference2DImg = reshape(differenceImg, imageSize);
        
        % Save to Library
        avgBkgLibrary(g).avgImageSubtracted = difference2DImg;

        %% Plot the Average 2D Image of Atom - Bkg for each Scan ID
        if analyVar.plotAvgAtomMinusBkg

            % Plot the 2D image
            figure;
            imagesc(difference2DImg);
            
            % Add a colorbar on the side
            colorbar;
            
            % Enhance the plot appearance (optional)
            colormap('turbo');       % Choose a colormap (e.g., 'jet', 'parula', 'hot')
            axis equal tight;      % Keep aspect ratio square and remove extra whitespace
            xlabel('X-axis');
            ylabel('Y-axis');
            title(sprintf( ...
            'Avg Image of Atom - Bkg: Scan ID %.4f and Scan Parameter %.4f', ...
            thisScanID,thisParam));
        end

    end

    %% 5. Save averaged images to a MAT file and TXT files

    outputDir = analyVar.analyOutDir;

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Stack averaged images as columns: pixels x groups.
    avgBkgImages = reshape( ...
        [avgBkgLibrary.avgImageBkg], numPixels, numGroups);
    avgAtomImages = reshape( ...
        [avgBkgLibrary.avgImageAtom], numPixels, numGroups);

    scanIDs = [avgBkgLibrary.scanID].';
    scanParameters = [avgBkgLibrary.scanParameter].';
    numImagesAveraged = [avgBkgLibrary.numImages].';

    % Save the complete library and its metadata in one MAT file.
    outputFiles.mat = fullfile(outputDir, 'avgBkgLibrary.mat');

    save(outputFiles.mat, ...
        'avgBkgLibrary', 'avgBkgImages', 'avgAtomImages', ...
        'scanIDs', 'scanParameters', 'numImagesAveraged', ...
        'imageSize', '-v7');

    % Save image data: each column is one averaged image.
    outputFiles.imagesTxtAtom = fullfile( ...
        outputDir, 'avgAtomLibrary_images.txt');

    writematrix(avgAtomImages, outputFiles.imagesTxtAtom, ...
        'Delimiter', '\t');

    outputFiles.imagesTxtBkg = fullfile( ...
        outputDir, 'avgBkgLibrary_images.txt');

    writematrix(avgBkgImages, outputFiles.imagesTxtBkg, ...
        'Delimiter', '\t');

    % Save group metadata separately.
    outputFiles.metadataTxt = fullfile( ...
        outputDir, 'avgBkgLibrary_metadata.txt');

    metadataTable = table( ...
        (1:numGroups).', scanIDs, scanParameters, ...
        numImagesAveraged, ...
        'VariableNames', ...
        {'GroupIndex', 'ScanID', 'ScanParameter', 'NumImages'});

    writetable(metadataTable, outputFiles.metadataTxt, ...
        'Delimiter', '\t');

    % Save the output files
    %avgBkgLibrary.outputFiles =  outputFiles;

end