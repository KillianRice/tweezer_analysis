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
basenameNum = 1;
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
fprintf('Using basenameNum = %d, image index k = %d\n', basenameNum, k);
fprintf('Default ROI radius = %g pixels\n\n', roiRadiusPix);

%% Ask user for ROI radius
userRadius = input(sprintf('Enter ROI radius in pixels [%g]: ', roiRadiusPix));

if ~isempty(userRadius)
    roiRadiusPix = userRadius;
end

%% Average all images in this file
avgPrelimRawAtoms = zeros(analyVar.matrixSize);

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
end

avgPrelimRawAtoms = avgPrelimRawAtoms ./ indivDataset{basenameNum}.CounterAtom;

%% Same ROI cut construction as check_pnts
roiSideLength = 2*analyVar.roiWinRadAtom(basenameNum) + 1;

[roiCutImage, roi_Index] = deal(zeros(roiSideLength));

roi_Index(:)   = indivDataset{basenameNum}.image_Index(indivDataset{basenameNum}.image_Index ~= 0);
roiCutImage(:) = prelimRawAtoms(indivDataset{basenameNum}.image_Index ~= 0);

Xpix = (-analyVar.roiWinRadAtom(basenameNum):analyVar.roiWinRadAtom(basenameNum)) ...
    + analyVar.cloudColCntrAtom(1);

Ypix = (-analyVar.roiWinRadAtom(basenameNum):analyVar.roiWinRadAtom(basenameNum)) ...
    + analyVar.cloudRowCntrAtom(1);

%% Click tweezer positions on ROI cut image
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
%%%Creating local coordinates within the ROI of the image
centersXY(:,1) = round(centersXY(:,1) - min(Xpix) + 1);
centersXY(:,2) = round(centersXY(:,2) - min(Ypix) + 1);

nTweezers = size(centersXY,1);

if nTweezers == 0
    warning('No tweezer centers selected.');
    tweezerROI = [];
    return;
end

fprintf('Selected %d tweezer square ROIs.\n', nTweezers);

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