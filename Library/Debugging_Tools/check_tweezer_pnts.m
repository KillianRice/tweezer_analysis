function tweezerROI = check_tweezer_pnts(varargin)
% Interactive tweezer ROI selector.
%
% Click tweezer centers on a representative image.
% Saves:
%   tweezerROI.centersXY  = [x y] pixel centers
%   tweezerROI.radiusPix  = ROI radius in pixels
%   tweezerROI.basenameNum
%   tweezerROI.imageIndex
%
% Output file:
%   analyVar.analyOutDir/tweezerROI.txt
%   analyVar.analyOutDir/tweezerROI.mat

%% Load variables
analyVar = AnalysisVariables;
indivDataset = get_indiv_batch_data(analyVar);

%% Defaults
%basenameNum = 1;
k = 1;
roiRadiusPix = 5;

if nargin >= 1
    basenameNum = varargin{1};
end

if nargin >= 2
    k = varargin{2};
end

if nargin >= 3
    roiRadiusPix = varargin{3};
end

fprintf('\nTweezer ROI selector\n');
%fprintf('Using basenameNum = %d, image index k = %d\n', basenameNum, k);
fprintf('Default ROI radius = %g pixels\n\n', roiRadiusPix);

%% Ask user for ROI radius
userRadius = input(sprintf('Enter ROI radius in pixels [%g]: ', roiRadiusPix));

if ~isempty(userRadius)
    roiRadiusPix = userRadius;
end

%% Average all image in the selected files
avgPrelimRawAtoms = zeros(analyVar.matrixSize);
num = 0;

for basenameNum = 1:analyVar.numBasenamesAtom
    for k = 1:indivDataset{basenameNum}.CounterAtom
    
        atomFile = [analyVar.dataDir char(indivDataset{basenameNum}.fileAtom(k)) analyVar.dataAtom];
        backFile = [analyVar.dataDir char(indivDataset{basenameNum}.fileBack(k)) analyVar.dataBack];
    
        sFID = fopen(atomFile,'rb','ieee-be');
        tFID = fopen(backFile,'rb','ieee-be');
    
        if sFID < 0
            error('Could not open atom file:\n%s', atomFile);
        end
    
        if tFID < 0
            error('Could not open background file:\n%s', backFile);
        end
    
        fullRawImageAtom = double(fread(sFID, analyVar.matrixSize, '*int16'));
        fullRawImageBack = double(fread(tFID, analyVar.matrixSize, '*int16'));
    
        fclose(sFID);
        fclose(tFID);
    
        if analyVar.UseImages_Fluorescence == 1
            prelimRawAtoms = fullRawImageAtom;
        else
            prelimRawAtoms = log(abs(fullRawImageBack)) - log(abs(fullRawImageAtom));
        end
    
        avgPrelimRawAtoms = avgPrelimRawAtoms + prelimRawAtoms;
        num = num + 1;
    end
end

avgPrelimRawAtoms = avgPrelimRawAtoms ./ num ;

%% Same ROI cut construction as check_pnts
roiSideLength = 2*analyVar.roiWinRadAtom(basenameNum) + 1;

[roiCutImage, roi_Index] = deal(zeros(roiSideLength));

roi_Index(:)   = indivDataset{basenameNum}.image_Index(indivDataset{basenameNum}.image_Index ~= 0);
roiCutImage(:) = avgPrelimRawAtoms(indivDataset{basenameNum}.image_Index ~= 0);

Xpix = (-analyVar.roiWinRadAtom(basenameNum):analyVar.roiWinRadAtom(basenameNum)) ...
    + analyVar.cloudColCntrAtom(1);

Ypix = (-analyVar.roiWinRadAtom(basenameNum):analyVar.roiWinRadAtom(basenameNum)) ...
    + analyVar.cloudRowCntrAtom(1);

%% Click tweezer positions on ROI cut image
%% Select tweezer positions

fprintf('\nHow would you like to define the tweezer ROIs?\n');
fprintf('  1 = Click centers on image\n');
fprintf('  2 = Paste ROI coordinates from tweezerROI.txt\n\n');

while true

    selection = input('Selection [1]: ');

    if isempty(selection)
        selection = 1;
    end

    if isscalar(selection) && ismember(selection,[1 2])
        break;
    end

    fprintf('Please enter 1 or 2.\n\n');

end


%% Mouse-click mode

if selection == 1

    figure(501); clf;

    pcolor(Ypix, Xpix, roiCutImage);
    shading flat;
    axis equal tight;
    colorbar;

    title({'Click tweezer centers on ROI cut image', ...
           'Press Enter when finished'});

    xlabel('Y pixel / row');
    ylabel('X pixel / column');

    fprintf('\nClick each tweezer center. Press Enter when finished.\n');

    [yClick, xClick] = ginput;

    centersXY = round([xClick(:), yClick(:)]);

    % Convert from image coordinates to local ROI coordinates
    centersXY(:,1) = round( ...
        centersXY(:,1) - min(Xpix) + 1);

    centersXY(:,2) = round( ...
        centersXY(:,2) - min(Ypix) + 1);

    % Use the ROI radius entered above
    nTweezers = size(centersXY,1);


%% Manual-entry mode

else

    fprintf('\nManual tweezer ROI entry\n');
    fprintf('Enter ROI data using the following format:\n\n');

    fprintf('1    51    44    5\n');
    fprintf('2    48    71    5\n');
    fprintf('3    27    67    5\n');
    fprintf('4    31    42    5\n\n');

    fprintf('Columns are:\n');
    fprintf('  index   xCenter   yCenter   roiHalfWidthPix\n\n');

    fprintf('Paste/type the complete set of rows below.\n');
    fprintf('Finish by entering a blank line.\n\n');

    %% Read rows from command window

    roiData = [];

    while true

        line = input('', 's');

        % Blank line = finished
        if isempty(strtrim(line))
            break;
        end

        % Convert this row to numbers
        rowData = sscanf(line,'%f').';

        % Validate row
        if numel(rowData) ~= 4

            fprintf('\n');
            fprintf('Invalid row:\n');
            fprintf('  %s\n\n',line);

            fprintf( ...
                'Expected: index  xCenter  yCenter  roiHalfWidthPix\n\n');

            continue;

        end

        roiData(end+1,:) = rowData; %#ok<AGROW>

    end

    %% Make sure something was entered

    if isempty(roiData)

        warning('No tweezer centers entered.');

        tweezerROI = [];

        return;

    end

    %% Extract coordinates

    centersXY = round(roiData(:,2:3));

    %% Extract ROI radius

    roiRadiusValues = round(roiData(:,4));

    nTweezers = size(centersXY,1);

    %% Validate indices

    expectedIndices = (1:nTweezers).';

    if any(round(roiData(:,1)) ~= expectedIndices)

        warning( ...
            ['ROI indices are not sequential. ' ...
             'They will be saved as 1 through %d.'], ...
            nTweezers);

    end

    %% Validate ROI radius

    if any(roiRadiusValues <= 0)

        error('ROI half-width values must be positive.');

    end

    %% Require same ROI radius

    if any(roiRadiusValues ~= roiRadiusValues(1))

        error( ...
            ['Different ROI half-widths were entered. ' ...
             'All ROIs must currently use the same half-width.']);

    end

    roiRadiusPix = roiRadiusValues(1);

    fprintf('\n');
    fprintf('Read %d tweezer ROIs.\n',nTweezers);

end


%% Check that ROIs were supplied

if nTweezers == 0

    warning('No tweezer centers selected.');

    tweezerROI = [];

    return;

end

fprintf( ...
    'Selected %d tweezer square ROIs.\n', ...
    nTweezers);


%% Check that ROIs were supplied

if nTweezers == 0

    warning('No tweezer centers selected.');

    tweezerROI = [];

    return;

end

fprintf( ...
    'Selected %d tweezer square ROIs.\n', ...
    nTweezers);

% figure(501); clf;
% 
% pcolor(Ypix, Xpix, roiCutImage);
% shading flat;
% axis equal tight;
% colorbar;
% 
% title({'Click tweezer centers on ROI cut image', ...
%        'Press Enter when finished'});
% 
% xlabel('Y pixel / row');
% ylabel('X pixel / column');
% 
% fprintf('\nClick each tweezer center. Press Enter when finished.\n');
% 
% [yClick, xClick] = ginput;
% 
% centersXY = round([xClick(:), yClick(:)]);
% %%%Creating local coordinates within the ROI of the image
% centersXY(:,1) = round(centersXY(:,1) - min(Xpix) + 1);
% centersXY(:,2) = round(centersXY(:,2) - min(Ypix) + 1);
% 
% nTweezers = size(centersXY,1);
% 
% if nTweezers == 0
%     warning('No tweezer centers selected.');
%     tweezerROI = [];
%     return;
% end
% 
% fprintf('Selected %d tweezer square ROIs.\n', nTweezers);
% 
%% Save ROI info
tweezerROI = struct;
tweezerROI.centersXY = centersXY;
tweezerROI.roiHalfWidthPix = roiRadiusPix;
tweezerROI.roiShape = 'square';
tweezerROI.basenameNum = basenameNum;
tweezerROI.imageIndex = k;
tweezerROI.atomFile = atomFile;
tweezerROI.backFile = backFile;

matFile = fullfile(analyVar.analyOutDir, 'tweezerROI.mat');
save(matFile, 'tweezerROI');

txtFile = fullfile(analyVar.analyOutDir, 'tweezerROI.txt');
fid = fopen(txtFile, 'w');

fprintf(fid, '# Tweezer ROI file\n');
fprintf(fid, '# ROI shape: square\n');
fprintf(fid, '# Columns: index xCenter yCenter roiHalfWidthPix\n');

for idx = 1:nTweezers
    fprintf(fid, '%d\t%d\t%d\t%d\n', ...
        idx, ...
        centersXY(idx,1), ...
        centersXY(idx,2), ...
        roiRadiusPix);
end

fclose(fid);

fprintf('Saved tweezer ROI data to:\n%s\n%s\n\n', matFile, txtFile);

%% Overlay preview on ROI cut image
figure(502); clf;

%pcolor(Ypix, Xpix, roiCutImage);
imagesc(roiCutImage)
axis image
shading flat;
%axis equal tight;
colorbar;
hold on;

title(sprintf('Selected Square Tweezer ROIs, half-width = %d px', roiRadiusPix));
xlabel('Y pixel / row');
ylabel('X pixel / column');

for idx = 1:nTweezers

    x0 = centersXY(idx,1);
    y0 = centersXY(idx,2);
    r  = roiRadiusPix;

    rectangle( ...
        'Position', [y0-r, x0-r, 2*r+1, 2*r+1], ...
        'EdgeColor', 'r', ...
        'LineWidth', 1.5);

    plot(y0, x0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);

    text(y0 + r + 1, x0, num2str(idx), ...
        'Color', 'w', ...
        'FontWeight', 'bold');
end

hold off;

%% Cropped ROI preview
nRows = ceil(sqrt(nTweezers));
nCols = ceil(nTweezers/nRows);

figure(503); clf;

for idx = 1:nTweezers
    
    x0 = centersXY(idx,1);
    y0 = centersXY(idx,2);
    r  = roiRadiusPix;
    
    roiImg = roiCutImage(x0-r:x0+r, ...
                        y0-r:y0+r);

    subplot(nRows, nCols, idx);
    imagesc(roiImg);
    axis image off;
    title(sprintf('ROI %d', idx));
end

sgtitle('Cropped Square Tweezer ROI Preview');

end