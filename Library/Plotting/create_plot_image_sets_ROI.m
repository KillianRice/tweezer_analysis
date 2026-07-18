function create_plot_image_sets_ROI(analyVar, indivDataset)
% Plot selected raw ROI images, full OD images, and tweezer OD cuts.
%
% Usage:
%   plot_image_sets_ROI
%   plot_image_sets_ROI(analyVar, indivDataset)
%
% The function prompts the user for:
%   - basename/file range
%   - image range within each file
%   - which image types to display
%
% Image types:
%   1. Raw atom ROI
%   2. Raw background ROI
%   3. Full saved OD ROI
%   4. Saved individual tweezer OD cuts


%% Load tweezer ROI coordinates
tweezerFile = fullfile(analyVar.analyOutDir, 'tweezerROI.mat');

if exist(tweezerFile,'file')
    roiData = load(tweezerFile, 'tweezerROI');
    tweezerROI = roiData.tweezerROI;
    showTweezerBoxes = true;
else
    tweezerROI = [];
    showTweezerBoxes = false;
    warning('No tweezerROI.mat found. Raw ROI plots will not show tweezer boxes.');
end

%% Determine available file range
numFiles = analyVar.numBasenamesAtom;

fprintf('\nROI image plotting utility\n');
fprintf('Available file indices: 1 through %d\n\n', numFiles);

%% Ask which files to plot
firstFile = input('First file index [1]: ');

if isempty(firstFile)
    firstFile = 1;
end

lastFile = input(sprintf('Last file index [%d]: ', firstFile));

if isempty(lastFile)
    lastFile = firstFile;
end

firstFile = round(firstFile);
lastFile = round(lastFile);

if firstFile < 1 || lastFile > numFiles || firstFile > lastFile
    error('Invalid file range: %d through %d.', firstFile, lastFile);
end

%% Ask which plot types to produce
fprintf('\nChoose plot types:\n');
fprintf('  1 = Raw atom ROI images\n');
fprintf('  2 = Raw background ROI images\n');
fprintf('  3 = Full saved OD ROI images\n');
fprintf('  4 = Individual tweezer OD cuts\n');
fprintf('You may enter several choices, for example: [1 3 4]\n\n');

plotChoices = input('Plot choices [[1 3 4]]: ');

if isempty(plotChoices)
    plotChoices = [1 3 4];
end

plotChoices = unique(plotChoices(:)');

validChoices = [1 2 3 4];

if any(~ismember(plotChoices, validChoices))
    error('Plot choices must contain only 1, 2, 3, or 4.');
end

%% Optional tweezer ROI information
plotTweezers = ismember(4, plotChoices);

if plotTweezers
    roiFile = fullfile(analyVar.analyOutDir, 'tweezerROI.mat');

    if ~exist(roiFile, 'file')
        error('Cannot find tweezer ROI file:\n%s', roiFile);
    end

    roiData = load(roiFile, 'tweezerROI');
    tweezerROI = roiData.tweezerROI;
    numTweezers = size(tweezerROI.centersXY,1);
else
    tweezerROI = [];
    numTweezers = 0;
end

%% Loop through selected files
for basenameNum = firstFile:lastFile

    numImages = indivDataset{basenameNum}.CounterAtom;

    fprintf('\nFile %d contains %d images.\n', basenameNum, numImages);

    firstImage = input('First image index [1]: ');

    if isempty(firstImage)
        firstImage = 1;
    end

    lastImage = input(sprintf('Last image index [%d]: ', ...
        min(numImages, firstImage + 8)));

    if isempty(lastImage)
        lastImage = min(numImages, firstImage + 8);
    end

    firstImage = round(firstImage);
    lastImage = round(lastImage);

    if firstImage < 1 || lastImage > numImages || firstImage > lastImage
        error('Invalid image range for file %d: %d through %d.', ...
            basenameNum, firstImage, lastImage);
    end

    selectedImages = firstImage:lastImage;
    numSelected = numel(selectedImages);

    %% Subplot layout
    nRows = ceil(sqrt(numSelected));
    nCols = ceil(numSelected/nRows);

    %% Plot Coordinates
    rowPix = (-analyVar.roiWinRadAtom(basenameNum): ...
           analyVar.roiWinRadAtom(basenameNum)) ...
           + analyVar.cloudRowCntrAtom(basenameNum);

    colPix = (-analyVar.roiWinRadAtom(basenameNum): ...
           analyVar.roiWinRadAtom(basenameNum)) ...
           + analyVar.cloudColCntrAtom(basenameNum);

    %% Plot raw atom ROI images
    if ismember(1, plotChoices)

        figAtom = figure;
        set(figAtom, 'Name', sprintf('Raw Atom ROI - File %d', basenameNum));

        atomAx = gobjects(numSelected,1);

        for plotIndex = 1:numSelected

            k = selectedImages(plotIndex);

            atomROI = reshape( ...
                 indivDataset{basenameNum}.AtomsCloud(:,k), ...
                 [1 1]*(2*analyVar.roiWinRadAtom(basenameNum)+1));

            atomAx(plotIndex) = subplot(nRows,nCols,plotIndex);
            imagesc(colPix, rowPix, atomROI);
            axis image;
            colorbar;
            hold on;
            
            xlabel('Column pixel');
            ylabel('Row pixel');
            
            overlay_tweezer_boxes_global(gca, tweezerROI, rowPix, colPix);
            
            
            hold off;

            title(make_image_title(indivDataset, basenameNum, k));
        end

        apply_common_color_limits(atomAx);
        sgtitle(sprintf('Raw Atom ROI Images - File %d', basenameNum));
    end

    %% Plot raw background ROI images
    if ismember(2, plotChoices)

        figBack = figure;
        set(figBack, 'Name', sprintf('Raw Background ROI - File %d', basenameNum));

        backAx = gobjects(numSelected,1);

        for plotIndex = 1:numSelected

            k = selectedImages(plotIndex);

            backROI = reshape( ...
                indivDataset{basenameNum}.BackgroundCloud(:,k), ...
                [1 1]*(2*analyVar.roiWinRadAtom(basenameNum)+1));

            backAx(plotIndex) = subplot(nRows,nCols,plotIndex);
            imagesc(colPix, rowPix, backROI);
            axis image;
            colorbar;
            hold on;
            
            xlabel('Column pixel');
            ylabel('Row pixel');
            
            overlay_tweezer_boxes_global(gca, tweezerROI, rowPix, colPix);
            
            hold off;

            title(make_image_title(indivDataset, basenameNum, k));
        end

        apply_common_color_limits(backAx);
        sgtitle(sprintf('Raw Background ROI Images - File %d', basenameNum));
    end

    %% Plot full saved OD ROI images
    if ismember(3, plotChoices)

        figOD = figure;
        set(figOD, 'Name', sprintf('Full OD ROI - File %d', basenameNum));

        odAx = gobjects(numSelected,1);

        for plotIndex = 1:numSelected

            k = selectedImages(plotIndex);

            odFile = [analyVar.analyOutDir ...
                char(indivDataset{basenameNum}.fileAtom(k)) ...
                analyVar.ODimageFilename];

            if ~exist(odFile,'file')
                warning('Missing full OD image:\n%s', odFile);

                fullOD = nan( ...
                    2*analyVar.roiWinRadAtom(basenameNum)+1);
            else
                fullOD = readmatrix(odFile, 'FileType','text');
            end

            odAx(plotIndex) = subplot(nRows,nCols,plotIndex);
            imagesc(colPix, rowPix, fullOD);
            axis image;
            colorbar;
            hold on;
            
            xlabel('Column pixel');
            ylabel('Row pixel');
            
            overlay_tweezer_boxes_global(gca, tweezerROI, rowPix, colPix);
            
            hold off;

            title(make_image_title(indivDataset, basenameNum, k));
        end

        apply_common_color_limits(odAx);
        sgtitle(sprintf('Full Saved OD ROI Images - File %d', basenameNum));
    end

    %% Plot saved individual tweezer OD cuts
    if ismember(4, plotChoices)

        % One figure per selected image so each figure contains all tweezers
        for plotIndex = 1:numSelected

            k = selectedImages(plotIndex);

            nTweezerRows = ceil(sqrt(numTweezers));
            nTweezerCols = ceil(numTweezers/nTweezerRows);

            figTweezer = figure;
            set(figTweezer, 'Name', sprintf( ...
                'Tweezer OD Cuts - File %d Image %d', basenameNum, k));

            tweezerAx = gobjects(numTweezers,1);

            for tweezerNum = 1:numTweezers

                odFile = [analyVar.analyOutDir ...
                    char(indivDataset{basenameNum}.fileAtom(k)) ...
                    sprintf('_Tweezer%03d', tweezerNum) ...
                    analyVar.ODimageFilename];

                if ~exist(odFile,'file')
                    warning('Missing tweezer OD image:\n%s', odFile);
                    tweezerOD = nan( ...
                        2*tweezerROI.roiHalfWidthPix+1);
                else
                    tweezerOD = readmatrix(odFile, 'FileType','text');
                end

                tweezerAx(tweezerNum) = subplot( ...
                    nTweezerRows, nTweezerCols, tweezerNum);

                imagesc(colPix, rowPix, tweezerOD);
                axis image off;
                colorbar;

                title(sprintf('Tweezer %d', tweezerNum));
            end

            apply_common_color_limits(tweezerAx);

            sgtitle(sprintf( ...
                'Tweezer OD Cuts - File %d, Image %d', ...
                basenameNum, k));
        end
    end
end

end



%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titleText = make_image_title(indivDataset, basenameNum, k)
% Create a useful subplot title.

if isfield(indivDataset{basenameNum}, 'imagevcoAtom') && ...
        numel(indivDataset{basenameNum}.imagevcoAtom) >= k

    titleText = sprintf('Image %d: %g', ...
        k, indivDataset{basenameNum}.imagevcoAtom(k));
else
    titleText = sprintf('Image %d', k);
end

end



function apply_common_color_limits(axHandles)
% Apply one common color scale to all valid axes in the current figure.

axHandles = axHandles(isgraphics(axHandles));

if isempty(axHandles)
    return;
end

climValues = get(axHandles, 'CLim');

if iscell(climValues)
    climValues = cell2mat(climValues);
end

commonCLim = [ ...
    min(climValues(:,1), [], 'omitnan'), ...
    max(climValues(:,2), [], 'omitnan')];

if all(isfinite(commonCLim)) && commonCLim(1) < commonCLim(2)
    set(axHandles, 'CLim', commonCLim);
end

end

function overlay_tweezer_boxes(axH, tweezerROI, imageSize)
% Overlay square tweezer ROIs on an ROI-local image.
%
% centersXY(:,1) is treated as row center.
% centersXY(:,2) is treated as column center.

if isempty(tweezerROI) || ~isfield(tweezerROI,'centersXY')
    return;
end

centersXY = tweezerROI.centersXY;
r = tweezerROI.roiHalfWidthPix;

hold(axH,'on');

for tweezerNum = 1:size(centersXY,1)

    rowCenter = centersXY(tweezerNum,1);
    colCenter = centersXY(tweezerNum,2);

    rowMin = rowCenter-r;
    rowMax = rowCenter+r;

    colMin = colCenter-r;
    colMax = colCenter+r;

    % Warn if ROI extends beyond the plotted image
    if rowMin < 1 || colMin < 1 || ...
            rowMax > imageSize(1) || colMax > imageSize(2)

        warning('Tweezer ROI %d extends beyond image bounds.', tweezerNum);
    end

    rectangle(axH, ...
        'Position', [colMin-0.5, rowMin-0.5, 2*r+1, 2*r+1], ...
        'EdgeColor', 'r', ...
        'LineWidth', 1.5);

    plot(axH, colCenter, rowCenter, ...
        'r+', ...
        'MarkerSize', 8, ...
        'LineWidth', 1.25);

    text(axH, colMax+1, rowCenter, num2str(tweezerNum), ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle');
end

hold(axH,'off');

end

function overlay_tweezer_boxes_global(axH, tweezerROI, rowPix, colPix)

centersXY = tweezerROI.centersXY;
r = tweezerROI.roiHalfWidthPix;

hold(axH,'on');

for tweezerNum = 1:size(centersXY,1)

    rowLocal = centersXY(tweezerNum,1);
    colLocal = centersXY(tweezerNum,2);

    if rowLocal < 1 || rowLocal > numel(rowPix) || ...
            colLocal < 1 || colLocal > numel(colPix)
        warning('Tweezer %d center is outside ROI axes.', tweezerNum);
        continue;
    end

    rowCenter = rowPix(rowLocal);
    colCenter = colPix(colLocal);

    rectangle(axH, ...
        'Position', [ ...
            colCenter-r-0.5, ...
            rowCenter-r-0.5, ...
            2*r+1, ...
            2*r+1], ...
        'EdgeColor','r', ...
        'LineWidth',1.5);

    plot(axH, colCenter, rowCenter, ...
        'r+', ...
        'MarkerSize',8, ...
        'LineWidth',1.25);

    text(axH, colCenter+r+1, rowCenter, num2str(tweezerNum), ...
        'Color','w', ...
        'FontWeight','bold');
end

hold(axH,'off');

end