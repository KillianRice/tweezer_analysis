% function [fit2DAxH, fit1DAxH, resAxH] = create_plot_fitEval_Averaged(analyVar, avgDataset)
% 
% %% Figure numbers
% figStruct.fig2DFit = analyVar.figNum.fig2DFit;
% figStruct.figRes   = analyVar.figNum.figRes;
% figStruct.fig1DFit = analyVar.figNum.fig1DFit;
% figStruct.fig1DBEC = analyVar.figNum.fig1DBEC;
% 
% figure(figStruct.fig2DFit);
% figure(figStruct.figRes);
% figure(figStruct.fig1DFit);
% 
% if strcmpi(analyVar.InitCase,'Bimodal')
%     figure(figStruct.fig1DBEC);
% end
% 
% %% Plot layout
% nPlots = avgDataset.CounterAtom;
% 
% if isfield(avgDataset,'SubPlotRows') && isfield(avgDataset,'SubPlotCols')
%     nRows = avgDataset.SubPlotRows;
%     nCols = avgDataset.SubPlotCols;
% else
%     nRows = ceil(sqrt(nPlots));
%     nCols = ceil(nPlots/nRows);
% end
% 
% %% Preallocate axes
% [fit2DAxH, fit1DAxH, resAxH] = deal(zeros(1,nPlots));
% 
% %% Choose reference basename index for window sizes
% basenameNum = 1;
% 
% for j = 1:nPlots
% 
%     %% Preallocate reconstructed fit image
%     roiDistImage = zeros(size(avgDataset.All_OD_Image{j}));
% 
%     %% Smooth averaged OD image
%     OD_Smooth = analyVar.smoothFilt(avgDataset.All_OD_Image{j}, analyVar.smoothFiltMat);
% 
%     %% Get ROI windows
%     if isfield(avgDataset,'roiWin_Index')
%         roiWin_Index = avgDataset.roiWin_Index;
%     else
%         error('avgDataset.roiWin_Index is missing. Save roiWin_Index into avgDataset during averaging.');
%     end
% 
%     %% Retrieve fitted OD windows
%     if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1
% 
%         OD_Image_Single = dlmread(avgDataset.avgODFiles{j});
%     
%         roiWin_Index = {true(size(OD_Image_Single))};
%         fitWinSize = size(OD_Image_Single,1);
% 
%         OD_Fit_ImageCell = cellfun( ...
%             @(x) reshape(OD_Image_Single(x), [1 1]*fitWinSize), ...
%             roiWin_Index, ...
%             'UniformOutput', 0);
%     
%     else
%         OD_Fit_ImageCell = cellfun( ...
%             @(x) reshape(avgDataset.All_OD_Image{j}(x), [1 1]*analyVar.funcFitWin(basenameNum)), ...
%             cellfun(analyVar.fitWinLogicInd, roiWin_Index, 'UniformOutput', 0), ...
%             'UniformOutput', 0);
%     end
% 
%     %% Retrieve fit coefficients
%     PCell = avgDataset.All_PCell{j};
% 
%     %% Weighting
%     errCell = cellfun(@(x,y) get_OD_weight(x(end), y), ...
%         PCell, OD_Fit_ImageCell, 'UniformOutput', 0);
% 
%     %% Generate fit model image
%     [Xgrid,Ygrid] = meshgrid(1:fitWinSize);%analyVar.funcFitWin(basenameNum));
% 
%     fitDistCell = cellfun( ...
%         @(x,y) reshape( ...
%             feval(str2func(analyVar.fitModel), ...
%             x(1:length(analyVar.InitCondList)), ...
%             [Xgrid(:), Ygrid(:), y(:)]), ...
%             [1 1]*fitWinSize), ... %analyVar.funcFitWin(basenameNum
%         PCell, errCell, 'UniformOutput', 0);
% 
%     %% Insert each fit window into full ROI image
%     for i = 1:sum(analyVar.LatticeAxesFit)
% 
%         roiDistTmp = roiWin_Index{i};
% 
%         roiDistTmp(analyVar.fitWinLogicInd(roiDistTmp)) = fitDistCell{i};
% 
%         roiDistImage = roiDistImage + roiDistTmp;
%     end
% 
%     %% Label for this averaged point
%     if isfield(avgDataset,'paramVals')
%         plotLabel = sprintf('%g %s', avgDataset.paramVals(j), analyVar.xDataUnit);
%     else
%         plotLabel = sprintf('%s %g',analyVar.avgScanParam, avgDataset.imagevcoAtom(j));
%     end
% 
%     %% Plot 2D fit
%     set(0,'CurrentFigure',figStruct.fig2DFit)
%     fit2DAxH(j) = subplot(nRows,nCols,j);
% 
%     pcolor(roiDistImage);
%     shading flat;
%     colorbar;
%     title(plotLabel);
%     grid off;
% 
%     if j == nPlots
%         set(gcf,'Name','2D Cloud Fit: Averaged Scans');
%         mtit('2D Cloud Fit: Averaged Scans','FontSize',16,'zoff',.025,'xoff',-.01);
%     end
% 
%     %% Plot residuals
%     set(0,'CurrentFigure',figStruct.figRes)
%     resAxH(j) = subplot(nRows,nCols,j);
% 
%     pcolor(avgDataset.All_OD_Image{j} - roiDistImage);
%     shading flat;
%     colorbar;
%     title(plotLabel);
%     grid off;
% 
%     if j == nPlots
%         set(gcf,'Name','Residuals: Averaged Scans');
%         mtit('Residuals: Averaged Scans','FontSize',16,'zoff',.025,'xoff',-.01);
%     end
% 
%     %% Plot 1D cross sections
%     set(0,'CurrentFigure',figStruct.fig1DFit)
% 
%     CloudCntr = round( ...
%         (analyVar.roiWinRadAtom(basenameNum) - ...
%         (analyVar.funcFitWin(basenameNum) - 1)/2) + ...
%         [avgDataset.All_fitParams{j}{1}.xCntr, ...
%          avgDataset.All_fitParams{j}{1}.yCntr]);
% 
%     % Clamp center to image bounds
%     CloudCntr(1) = max(1, min(size(avgDataset.All_OD_Image{j},2), CloudCntr(1)));
%     CloudCntr(2) = max(1, min(size(avgDataset.All_OD_Image{j},1), CloudCntr(2)));
% 
%     OD_1D_Xdata = OD_Smooth(:,CloudCntr(1));
%     OD_1D_Ydata = OD_Smooth(CloudCntr(2),:);
% 
%     OD_1D_Xfit = roiDistImage(:,CloudCntr(1));
%     OD_1D_Yfit = roiDistImage(CloudCntr(2),:);
% 
%     fit1DAxH(j) = subplot(nRows,nCols,j);
% 
%     hold on;
%     grid off;
% 
%     plot(OD_1D_Xdata,'c.');
%     plot(OD_1D_Ydata,'g.');
%     plot(OD_1D_Xfit,'k');
%     plot(OD_1D_Yfit,'r');
% 
%     title(sprintf('%s [%g,%g]', plotLabel, CloudCntr(1), CloudCntr(2)));
% 
%     xlim(mean(CloudCntr) + [-1 1]*analyVar.roiWinRadAtom(basenameNum));
% 
%     ylim([-0.1, ...
%         max(max([OD_1D_Xdata; OD_1D_Ydata'; OD_1D_Xfit; OD_1D_Yfit'])) + 0.1]);
% 
%     if j == nPlots
%         legend('y data','x data','y fit','x fit');
%         set(gcf,'Name','1D Fit: Averaged Scans');
% 
%         mtit(['Cross-Section of Averaged Fit using ' ...
%             strrep(analyVar.fitModel, analyVar.InitCase, [analyVar.InitCase ' '])], ...
%             'FontSize',16,'zoff',.025,'xoff',-.01);
%     end
% 
% end
% 
% end

function [fit2DAxH,fit1DAxH,resAxH] = ...
    create_plot_fitEval_Averaged(analyVar,avgDataset)
% Plot fitted averaged OD images.
%
% Supports:
%
%   avgDataset.AveragedFitMode = 'scanParameter'
%
%   avgDataset.AveragedFitMode = 'allImagesByTweezer'

%% Determine plotting mode
if isfield(avgDataset,'AveragedFitMode') && ...
        ~isempty(avgDataset.AveragedFitMode)

    fitMode = avgDataset.AveragedFitMode;

elseif isfield(analyVar,'AveragedFitMode') && ...
        ~isempty(analyVar.AveragedFitMode)

    fitMode = analyVar.AveragedFitMode;

else

    fitMode = 'scanParameter';

end

%% Build a common plotting collection
plotData = get_averaged_fit_plot_data( ...
    analyVar,avgDataset,fitMode);

nPlots = plotData.nPlots;

if nPlots == 0
    error('No fitted averaged images are available to plot.');
end

%% Figure numbers
figStruct.fig2DFit = analyVar.figNum.fig2DFit;
figStruct.figRes = analyVar.figNum.figRes;
figStruct.fig1DFit = analyVar.figNum.fig1DFit;
figStruct.fig1DBEC = analyVar.figNum.fig1DBEC;
figStruct.avgODImages = analyVar.figNum.avgODImages;

figure(figStruct.avgODImages);
clf;

figure(figStruct.fig2DFit);
clf;

figure(figStruct.figRes);
clf;

figure(figStruct.fig1DFit);
clf;

if strcmpi(analyVar.InitCase,'Bimodal')
    figure(figStruct.fig1DBEC);
end

%% Plot layout
if strcmpi(fitMode,'scanParameter') && ...
        isfield(avgDataset,'SubPlotRows') && ...
        isfield(avgDataset,'SubPlotCols') && ...
        avgDataset.SubPlotRows*avgDataset.SubPlotCols >= nPlots

    nRows = avgDataset.SubPlotRows;
    nCols = avgDataset.SubPlotCols;

else

    nRows = floor(sqrt(nPlots));
    nRows = max(nRows,1);

    nCols = ceil(nPlots/nRows);

end

%% Preallocate axes
avgImageAxH = gobjects(1,nPlots);
fit2DAxH = gobjects(1,nPlots);
fit1DAxH = gobjects(1,nPlots);
resAxH = gobjects(1,nPlots);

%% Store images for global color scales
fitImageCell = cell(nPlots,1);
residualImageCell = cell(nPlots,1);

%% First pass: reconstruct all fitted images
for j = 1:nPlots

    OD_Image = plotData.images{j};
    PCell = plotData.PCell{j};

    if isempty(OD_Image) || isempty(PCell)
        continue;
    end

    basenameNum = plotData.referenceBasenameNums(j);

    if isnan(basenameNum) || basenameNum < 1
        basenameNum = 1;
    end

    [roiDistImage,OD_Fit_ImageCell] = ...
        reconstruct_fit_image( ...
            analyVar, ...
            OD_Image, ...
            PCell, ...
            basenameNum);

    fitImageCell{j} = roiDistImage;
    residualImageCell{j} = OD_Image-roiDistImage;

    plotData.fitWindowCells{j} = OD_Fit_ImageCell;

end

%% Compute common color limits
fitCLim = get_global_image_limits(fitImageCell);
resCLim = get_symmetric_global_image_limits(residualImageCell);
avgImageCLim = get_global_image_limits(plotData.images);

%% Plot each result
for j = 1:nPlots

    OD_Image = plotData.images{j};
    PCell = plotData.PCell{j};
    roiDistImage = fitImageCell{j};

    if isempty(OD_Image) || ...
            isempty(PCell) || ...
            isempty(roiDistImage)

        continue;

    end

    residualImage = residualImageCell{j};

    basenameNum = plotData.referenceBasenameNums(j);

    if isnan(basenameNum) || basenameNum < 1
        basenameNum = 1;
    end

    plotLabel = plotData.labels{j};

    %% Plot saved averaged OD image
    set(0,'CurrentFigure',figStruct.avgODImages);

    avgImageAxH(j) = subplot(nRows,nCols,j);

    imagesc(OD_Image);
    set(avgImageAxH(j),'YDir','normal');

    axis(avgImageAxH(j),'image');
    axis(avgImageAxH(j),'tight');

    colorbar(avgImageAxH(j));

    % Apply the same color scale to every averaged OD image
    if all(isfinite(avgImageCLim)) && ...
            avgImageCLim(1) < avgImageCLim(2)

        clim(avgImageAxH(j),avgImageCLim);

    end

    xlabel(avgImageAxH(j),'Column pixel');
    ylabel(avgImageAxH(j),'Row pixel');

    title( ...
        avgImageAxH(j), ...
        plotLabel, ...
        'Interpreter','none');

    grid(avgImageAxH(j),'off');
    box(avgImageAxH(j),'on');

    %% Smoothed image for 1D display
    if analyVar.fitSmoothOD

        OD_Smooth = analyVar.smoothFilt( ...
            OD_Image, ...
            analyVar.smoothFiltMat);

    else

        OD_Smooth = OD_Image;

    end

    %% Plot 2D fitted distribution
    set(0,'CurrentFigure',figStruct.fig2DFit);

    fit2DAxH(j) = subplot(nRows,nCols,j);

    imagesc(roiDistImage);
    set(gca,'YDir','normal');

    axis image tight;

    colorbar;

    if all(isfinite(fitCLim)) && fitCLim(1) < fitCLim(2)
        clim(fit2DAxH(j),fitCLim);
    end

    xlabel('Column pixel');
    ylabel('Row pixel');

    title(plotLabel,'Interpreter','none');

    grid off;
    box on;

    %% Plot residual
    set(0,'CurrentFigure',figStruct.figRes);

    resAxH(j) = subplot(nRows,nCols,j);

    imagesc(residualImage);
    set(gca,'YDir','normal');

    axis image tight;

    colorbar;

    if all(isfinite(resCLim)) && resCLim(1) < resCLim(2)
        clim(resAxH(j),resCLim);
    end

    xlabel('Column pixel');
    ylabel('Row pixel');

    title(plotLabel,'Interpreter','none');

    grid off;
    box on;

    %% Determine fitted cloud center
    CloudCntr = get_fit_center_from_pcell( ...
        analyVar, ...
        PCell, ...
        size(OD_Image), ...
        basenameNum, ...
        strcmpi(fitMode,'allImagesByTweezer') || ...
            (isfield(analyVar,'UseTweezer') && ...
             analyVar.UseTweezer == 1));

    %% Clamp center to image bounds
    CloudCntr(1) = max( ...
        1,min(size(OD_Image,2),round(CloudCntr(1))));

    CloudCntr(2) = max( ...
        1,min(size(OD_Image,1),round(CloudCntr(2))));

    %% Extract 1D cross sections
    OD_1D_Xdata = OD_Smooth(:,CloudCntr(1));
    OD_1D_Ydata = OD_Smooth(CloudCntr(2),:);

    OD_1D_Xfit = roiDistImage(:,CloudCntr(1));
    OD_1D_Yfit = roiDistImage(CloudCntr(2),:);

    %% Plot 1D cross sections
    set(0,'CurrentFigure',figStruct.fig1DFit);

    fit1DAxH(j) = subplot(nRows,nCols,j);

    hold on;

    plot(OD_1D_Xdata,'c.');
    plot(OD_1D_Ydata,'g.');

    plot(OD_1D_Xfit,'k');
    plot(OD_1D_Yfit,'r');

    title( ...
        sprintf( ...
            '%s, center [%d,%d]', ...
            plotLabel, ...
            CloudCntr(1), ...
            CloudCntr(2)), ...
        'Interpreter','none');

    xlabel('Pixel');
    ylabel('OD');

    grid off;
    box on;

    yValues = [
        OD_1D_Xdata(:)
        OD_1D_Ydata(:)
        OD_1D_Xfit(:)
        OD_1D_Yfit(:)
        ];

    validY = yValues(isfinite(yValues));

    if ~isempty(validY)

        yPadding = 0.05*max( ...
            max(validY)-min(validY), ...
            1e-6);

        ylim([
            min(validY)-yPadding
            max(validY)+yPadding
            ]);

    end

    hold off;

end

%% Figure titles
set(0,'CurrentFigure',figStruct.avgODImages);

if strcmpi(fitMode,'allImagesByTweezer')

    set( ...
        gcf, ...
        'Name', ...
        'Averaged OD Images: All Individual Images by Tweezer');

    mtit( ...
        'Averaged OD Images: All Individual Images by Tweezer', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

else

    set( ...
        gcf, ...
        'Name', ...
        'Averaged OD Images: Scan Parameters');

    mtit( ...
        'Averaged OD Images: Scan Parameters', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

end

set(0,'CurrentFigure',figStruct.fig2DFit);

if strcmpi(fitMode,'allImagesByTweezer')

    set(gcf,'Name','2D Cloud Fit: All Images Averaged by Tweezer');

    mtit( ...
        '2D Cloud Fit: All Images Averaged by Tweezer', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

else

    set(gcf,'Name','2D Cloud Fit: Averaged Scans');

    mtit( ...
        '2D Cloud Fit: Averaged Scans', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

end

set(0,'CurrentFigure',figStruct.figRes);

if strcmpi(fitMode,'allImagesByTweezer')

    set(gcf,'Name','Residuals: All Images Averaged by Tweezer');

    mtit( ...
        'Residuals: All Images Averaged by Tweezer', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

else

    set(gcf,'Name','Residuals: Averaged Scans');

    mtit( ...
        'Residuals: Averaged Scans', ...
        'FontSize',16, ...
        'zoff',.025, ...
        'xoff',-.01);

end

set(0,'CurrentFigure',figStruct.fig1DFit);

if strcmpi(fitMode,'allImagesByTweezer')

    oneDTitle = [
        'Cross-Sections: All Images Averaged by Tweezer using ' ...
        strrep( ...
            analyVar.fitModel, ...
            analyVar.InitCase, ...
            [analyVar.InitCase ' '])
        ];

    set(gcf,'Name','1D Fit: All Images Averaged by Tweezer');

else

    oneDTitle = [
        'Cross-Section of Averaged Fit using ' ...
        strrep( ...
            analyVar.fitModel, ...
            analyVar.InitCase, ...
            [analyVar.InitCase ' '])
        ];

    set(gcf,'Name','1D Fit: Averaged Scans');

end

mtit( ...
    oneDTitle, ...
    'FontSize',16, ...
    'zoff',.025, ...
    'xoff',-.01);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData = get_averaged_fit_plot_data( ...
    analyVar,avgDataset,fitMode)

plotData = struct;

switch lower(fitMode)

    case lower('allImagesByTweezer')

        if ~isfield(avgDataset,'tweezerAllImageFit')

            error([ ...
                'avgDataset.tweezerAllImageFit is missing.\n' ...
                'Run imagefit_NumDistFit_Averaged with:\n' ...
                'analyVar.AveragedFitMode = ''allImagesByTweezer'';']);

        end

        fitData = avgDataset.tweezerAllImageFit;

        plotData.images = fitData.images;
        plotData.PCell = fitData.PCell;

        plotData.nPlots = fitData.numTweezers;

        plotData.referenceBasenameNums = ...
            fitData.referenceBasenameNums;

        plotData.referenceImageNums = ...
            fitData.referenceImageNums;

        plotData.labels = cell(plotData.nPlots,1);

        for tweezerNum = 1:plotData.nPlots

            numImages = ...
                fitData.numImagesAveraged(tweezerNum);

            plotData.labels{tweezerNum} = sprintf( ...
                'Tweezer %d, N = %d', ...
                tweezerNum,numImages);

        end

    otherwise

        if ~isfield(avgDataset,'All_OD_Image')
            error('avgDataset.All_OD_Image is missing.');
        end

        if ~isfield(avgDataset,'All_PCell')
            error('avgDataset.All_PCell is missing.');
        end

        plotData.images = avgDataset.All_OD_Image;
        plotData.PCell = avgDataset.All_PCell;

        plotData.nPlots = min( ...
            numel(plotData.images), ...
            numel(plotData.PCell));

        plotData.referenceBasenameNums = ...
            ones(plotData.nPlots,1);

        plotData.referenceImageNums = ...
            ones(plotData.nPlots,1);

        plotData.labels = cell(plotData.nPlots,1);

        for j = 1:plotData.nPlots

            if isfield(avgDataset,'sourceInfo') && ...
                    numel(avgDataset.sourceInfo) >= j && ...
                    ~isempty(avgDataset.sourceInfo{j})

                src = avgDataset.sourceInfo{j};

                plotData.referenceBasenameNums(j) = ...
                    src.basenameNums(1);

                plotData.referenceImageNums(j) = ...
                    src.imageNums(1);

                if isfield(src,'tweezerNums') && ...
                        ~isempty(src.tweezerNums)

                    tweezerNum = src.tweezerNums(1);

                else

                    tweezerNum = get_tweezer_number_for_plot( ...
                        j,avgDataset);

                end

            else

                tweezerNum = get_tweezer_number_for_plot( ...
                    j,avgDataset);

            end

            if isfield(avgDataset,'paramVals') && ...
                    numel(avgDataset.paramVals) >= j

                parameterValue = avgDataset.paramVals(j);

            elseif isfield(avgDataset,'imagevcoAtom') && ...
                    numel(avgDataset.imagevcoAtom) >= j

                parameterValue = avgDataset.imagevcoAtom(j);

            else

                parameterValue = j;

            end

            plotData.labels{j} = sprintf( ...
                'Tweezer %d, %s = %g %s', ...
                tweezerNum, ...
                analyVar.avgScanParam, ...
                parameterValue, ...
                analyVar.xDataUnit);

        end

end

plotData.fitWindowCells = cell(plotData.nPlots,1);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [roiDistImage,OD_Fit_ImageCell] = ...
    reconstruct_fit_image( ...
        analyVar,OD_Image,PCell,basenameNum)

OD_Image = double(OD_Image);

%% Optional smoothing is used only as model weighting input
if analyVar.fitSmoothOD

    OD_FitInput = analyVar.smoothFilt( ...
        OD_Image, ...
        analyVar.smoothFiltMat);

else

    OD_FitInput = OD_Image;

end

%% Tweezer image is already a complete fit window
if isfield(analyVar,'UseTweezer') && ...
        analyVar.UseTweezer == 1

    fitWinSize = size(OD_Image,1);

    if size(OD_Image,2) ~= fitWinSize
        error('Tweezer OD image must be square.');
    end

    OD_Fit_ImageCell = {
        OD_FitInput
        };

    roiWin_Index = {
        true(size(OD_Image))
        };

else

    error([ ...
        'The replacement plotting helper currently expects saved ' ...
        'individual tweezer cuts when UseTweezer is enabled.\n' ...
        'The original non-tweezer reconstruction path can be retained ' ...
        'if full-cloud fitting is still required.']);

end

%% Weighting
errCell = cellfun( ...
    @(fitParameters,fitImage) ...
        get_OD_weight( ...
            fitParameters(end),fitImage), ...
    PCell, ...
    OD_Fit_ImageCell, ...
    'UniformOutput',false);

%% Model grid
[Xgrid,Ygrid] = meshgrid(1:fitWinSize);

%% Evaluate model
fitDistCell = cellfun( ...
    @(fitParameters,fitError) reshape( ...
        feval( ...
            str2func(analyVar.fitModel), ...
            fitParameters( ...
                1:length(analyVar.InitCondList)), ...
            [ ...
                Xgrid(:), ...
                Ygrid(:), ...
                fitError(:)]), ...
        [fitWinSize fitWinSize]), ...
    PCell, ...
    errCell, ...
    'UniformOutput',false);

%% Insert fit into image
if numel(fitDistCell) == 1 && ...
        all(roiWin_Index{1}(:))

    roiDistImage = fitDistCell{1};

else

    roiDistImage = zeros(size(OD_Image));

    for fitNum = 1:numel(fitDistCell)

        roiMask = roiWin_Index{fitNum};

        roiMask(roiMask) = fitDistCell{fitNum};

        roiDistImage = roiDistImage+roiMask;

    end

end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloudCntr = get_fit_center_from_pcell( ...
    analyVar,PCell,imageSize,basenameNum,isLocalTweezerImage)

fitParameters = PCell{1};

parameterNames = analyVar.InitCondList;

if isstring(parameterNames)
    parameterNames = cellstr(parameterNames);
end

xIndex = find( ...
    strcmpi(parameterNames,'xCntr') | ...
    strcmpi(parameterNames,'xCenter') | ...
    strcmpi(parameterNames,'x0'), ...
    1);

yIndex = find( ...
    strcmpi(parameterNames,'yCntr') | ...
    strcmpi(parameterNames,'yCenter') | ...
    strcmpi(parameterNames,'y0'), ...
    1);

if isempty(xIndex) || isempty(yIndex)

    CloudCntr = [
        (imageSize(2)+1)/2
        (imageSize(1)+1)/2
        ];

    return;

end

fitX = fitParameters(xIndex);
fitY = fitParameters(yIndex);

if isLocalTweezerImage

    % Fit coordinates already correspond to the saved tweezer cut.
    CloudCntr = [fitX fitY];

else

    % Preserve the original full-ROI offset convention.
    CloudCntr = ...
        (analyVar.roiWinRadAtom(basenameNum) - ...
        (analyVar.funcFitWin(basenameNum)-1)/2) + ...
        [fitX fitY];

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%
function colorLimits = get_global_image_limits(imageCell)

globalMin = inf;
globalMax = -inf;

for j = 1:numel(imageCell)

    imageData = imageCell{j};

    if isempty(imageData)
        continue;
    end

    validData = imageData(isfinite(imageData));

    if isempty(validData)
        continue;
    end

    globalMin = min(globalMin,min(validData));
    globalMax = max(globalMax,max(validData));

end

if isfinite(globalMin) && ...
        isfinite(globalMax) && ...
        globalMin < globalMax

    colorLimits = [globalMin globalMax];

else

    colorLimits = [NaN NaN];

end

end

function colorLimits = ...
    get_symmetric_global_image_limits(imageCell)

largestMagnitude = 0;

for j = 1:numel(imageCell)

    imageData = imageCell{j};

    if isempty(imageData)
        continue;
    end

    validData = imageData(isfinite(imageData));

    if isempty(validData)
        continue;
    end

    largestMagnitude = max( ...
        largestMagnitude, ...
        max(abs(validData)));

end

if largestMagnitude > 0

    colorLimits = [
        -largestMagnitude
         largestMagnitude
        ];

else

    colorLimits = [NaN NaN];

end

end


function tweezerNum = ...
    get_tweezer_number_for_plot(j,avgDataset)

if isfield(avgDataset,'numTweezers') && ...
        avgDataset.numTweezers > 0

    tweezerNum = ...
        mod(j-1,avgDataset.numTweezers)+1;

else

    tweezerNum = 1;

end

end