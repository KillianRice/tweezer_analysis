function [fit2DAxH,fit1DAxH,resAxH] = ...
    create_plot_fitEval_Averaged(analyVar,avgDataset)
% Plot fitted averaged OD images.
%
% This function plots the averaged images produced by
% imagefit_NumDistFit_Averaged in scan-parameter mode.
%
% For each averaged image, it creates:
%   1. Averaged OD image
%   2. 2D fitted distribution
%   3. Residual image
%   4. 1D X/Y cross-sections of data and fit
%
% Required avgDataset fields:
%   avgDataset.All_OD_Image
%   avgDataset.All_PCell
%
% Optional fields used for labeling/layout:
%   avgDataset.sourceInfo
%   avgDataset.paramVals
%   avgDataset.imagevcoAtom
%   avgDataset.SubPlotRows
%   avgDataset.SubPlotCols


%% ================================================================
%  1. COLLECT DATA TO BE PLOTTED
%  ================================================================

plotData = get_averaged_fit_plot_data( ...
    analyVar, ...
    avgDataset);

nPlots = plotData.nPlots;

if nPlots == 0

    error( ...
        'No fitted averaged images are available to plot.');

end


%% ================================================================
%  2. FIGURE NUMBERS
%  ================================================================

figStruct.fig2DFit   = analyVar.figNum.fig2DFit;
figStruct.figRes      = analyVar.figNum.figRes;
figStruct.fig1DFit    = analyVar.figNum.fig1DFit;
figStruct.fig1DBEC    = analyVar.figNum.fig1DBEC;
figStruct.avgODImages = analyVar.figNum.avgODImages;


% Clear/reuse the existing figures.

figure(figStruct.avgODImages);
clf;

figure(figStruct.fig2DFit);
clf;

figure(figStruct.figRes);
clf;

figure(figStruct.fig1DFit);
clf;


% Only create/use the BEC figure for bimodal fits.

if strcmpi(analyVar.InitCase,'Bimodal')

    figure(figStruct.fig1DBEC);

end


%% ================================================================
%  3. DETERMINE SUBPLOT LAYOUT
%  ================================================================

if isfield(avgDataset,'SubPlotRows') && ...
        isfield(avgDataset,'SubPlotCols') && ...
        avgDataset.SubPlotRows * avgDataset.SubPlotCols >= nPlots

    nRows = avgDataset.SubPlotRows;
    nCols = avgDataset.SubPlotCols;

else

    % Automatically choose a roughly square layout.

    nRows = floor(sqrt(nPlots));
    nRows = max(nRows,1);

    nCols = ceil(nPlots/nRows);

end


%% ================================================================
%  4. PREALLOCATE AXES HANDLES
%  ================================================================

avgImageAxH = gobjects(1,nPlots);
fit2DAxH    = gobjects(1,nPlots);
fit1DAxH    = gobjects(1,nPlots);
resAxH      = gobjects(1,nPlots);


%% ================================================================
%  5. PREPARE STORAGE FOR GLOBAL COLOR LIMITS
%  ================================================================
%
% The color scale will be the same for every image in a given
% figure. This makes it easier to compare the different scan
% parameter values.

fitImageCell      = cell(nPlots,1);
residualImageCell = cell(nPlots,1);


%% ================================================================
%  6. FIRST PASS: RECONSTRUCT ALL FITTED IMAGES
%  ================================================================
%
% We do this before plotting so that the global color limits can
% be calculated from every image.

for j = 1:nPlots

    OD_Image = plotData.images{j};
    PCell    = plotData.PCell{j};


    % Skip entries that do not contain valid data.

    if isempty(OD_Image) || isempty(PCell)

        continue;

    end


    %% Determine which original image geometry to use

    basenameNum = ...
        plotData.referenceBasenameNums(j);

    if isnan(basenameNum) || basenameNum < 1

        basenameNum = 1;

    end


    %% Reconstruct the fitted 2D distribution

    [roiDistImage,OD_Fit_ImageCell] = ...
        reconstruct_fit_image( ...
            analyVar, ...
            OD_Image, ...
            PCell, ...
            basenameNum);


    %% Store fit and residual

    fitImageCell{j} = roiDistImage;

    residualImageCell{j} = ...
        OD_Image - roiDistImage;


    %% Keep the fit-window information

    plotData.fitWindowCells{j} = ...
        OD_Fit_ImageCell;

end


%% ================================================================
%  7. CALCULATE GLOBAL COLOR LIMITS
%  ================================================================

% Color scale for fitted images.

fitCLim = ...
    get_global_image_limits(fitImageCell);


% Symmetric color scale for residuals.

resCLim = ...
    get_symmetric_global_image_limits( ...
        residualImageCell);


% Color scale for averaged OD images.

avgImageCLim = ...
    get_global_image_limits( ...
        plotData.images);


%% ================================================================
%  8. PLOT EACH AVERAGED RESULT
%  ================================================================

for j = 1:nPlots

    %% ------------------------------------------------------------
    % Retrieve data for this scan point
    % -------------------------------------------------------------

    OD_Image = plotData.images{j};
    PCell    = plotData.PCell{j};

    roiDistImage = fitImageCell{j};


    % Empty data should normally have been caught during the first
    % pass. Throw an explicit error here so the problem is obvious.

    if isempty(OD_Image) || isempty(PCell)

        error( ...
            'Plot %d has empty data: OD_Image=%d, PCell=%d', ...
            j, ...
            isempty(OD_Image), ...
            isempty(PCell));

    end


    residualImage = ...
        residualImageCell{j};


    %% ------------------------------------------------------------
    % Recover original-image geometry
    % -------------------------------------------------------------

    basenameNum = ...
        plotData.referenceBasenameNums(j);

    if isnan(basenameNum) || basenameNum < 1

        basenameNum = 1;

    end


    %% ------------------------------------------------------------
    % Get plot label
    % -------------------------------------------------------------

    plotLabel = ...
        plotData.labels{j};


    %% ============================================================
    % 8A. PLOT AVERAGED OD IMAGE
    % =============================================================

    set( ...
        0, ...
        'CurrentFigure', ...
        figStruct.avgODImages);


    avgImageAxH(j) = ...
        subplot(nRows,nCols,j);


    imagesc(OD_Image);

    set( ...
        avgImageAxH(j), ...
        'YDir', ...
        'normal');


    axis( ...
        avgImageAxH(j), ...
        'image');

    axis( ...
        avgImageAxH(j), ...
        'tight');


    colorbar( ...
        avgImageAxH(j));


    % Apply common color scale.

    if all(isfinite(avgImageCLim)) && ...
            avgImageCLim(1) < avgImageCLim(2)

        clim( ...
            avgImageAxH(j), ...
            avgImageCLim);

    end


    xlabel( ...
        avgImageAxH(j), ...
        'Column pixel');

    ylabel( ...
        avgImageAxH(j), ...
        'Row pixel');


    title( ...
        avgImageAxH(j), ...
        plotLabel, ...
        'Interpreter','none');


    grid( ...
        avgImageAxH(j), ...
        'off');

    box( ...
        avgImageAxH(j), ...
        'on');


    %% ============================================================
    % 8B. CREATE SMOOTHED IMAGE FOR 1D DISPLAY
    % =============================================================

    if analyVar.fitSmoothOD

        OD_Smooth = ...
            analyVar.smoothFilt( ...
                OD_Image, ...
                analyVar.smoothFiltMat);

    else

        OD_Smooth = OD_Image;

    end


    %% ============================================================
    % 8C. PLOT 2D FIT
    % =============================================================

    set( ...
        0, ...
        'CurrentFigure', ...
        figStruct.fig2DFit);


    fit2DAxH(j) = ...
        subplot(nRows,nCols,j);


    imagesc(roiDistImage);

    set( ...
        fit2DAxH(j), ...
        'YDir', ...
        'normal');


    axis image tight;

    colorbar;


    % Apply common fit color scale.

    if all(isfinite(fitCLim)) && ...
            fitCLim(1) < fitCLim(2)

        clim( ...
            fit2DAxH(j), ...
            fitCLim);

    end


    xlabel('Column pixel');
    ylabel('Row pixel');

    title( ...
        plotLabel, ...
        'Interpreter','none');

    grid off;
    box on;


    %% ============================================================
    % 8D. PLOT RESIDUAL
    % =============================================================

    set( ...
        0, ...
        'CurrentFigure', ...
        figStruct.figRes);


    resAxH(j) = ...
        subplot(nRows,nCols,j);


    imagesc(residualImage);

    set( ...
        resAxH(j), ...
        'YDir', ...
        'normal');


    axis image tight;

    colorbar;


    % Use one symmetric color scale for all residuals.

    if all(isfinite(resCLim)) && ...
            resCLim(1) < resCLim(2)

        clim( ...
            resAxH(j), ...
            resCLim);

    end


    xlabel('Column pixel');
    ylabel('Row pixel');

    title( ...
        plotLabel, ...
        'Interpreter','none');

    grid off;
    box on;


    %% ============================================================
    % 8E. DETERMINE FITTED CLOUD CENTER
    % =============================================================

    CloudCntr = ...
        get_fit_center_from_pcell( ...
            analyVar, ...
            PCell, ...
            size(OD_Image), ...
            basenameNum, ...
            isfield(analyVar,'UseTweezer') && ...
            analyVar.UseTweezer == 1);


    %% Clamp center to image boundaries

    CloudCntr(1) = max( ...
        1, ...
        min( ...
            size(OD_Image,2), ...
            round(CloudCntr(1))));


    CloudCntr(2) = max( ...
        1, ...
        min( ...
            size(OD_Image,1), ...
            round(CloudCntr(2))));


    %% ============================================================
    % 8F. EXTRACT 1D DATA AND FIT CROSS-SECTIONS
    % =============================================================

    % X direction / column through cloud center.
    
    if isfield(analyVar,'dimenReduc1D') && ...
            analyVar.dimenReduc1D == 1
    
        %% ------------------------------------------------------------
        %  Integrated profiles
        % ------------------------------------------------------------
    
        OD_1D_Xdata = sum(OD_Smooth,1);
        OD_1D_Ydata = sum(OD_Smooth,2)';
    
        %% Evaluate the fitted 2D model and integrate it
    
        OD_1D_Xfit = sum(roiDistImage,1);
        OD_1D_Yfit = sum(roiDistImage,2)';
    
    else
    
        %% ------------------------------------------------------------
        %  Existing center cross-sections
        % ------------------------------------------------------------
    
        OD_1D_Xdata = ...
            OD_Smooth(:,CloudCntr(1));
    
        OD_1D_Xfit = ...
            roiDistImage(:,CloudCntr(1));
    
    
        OD_1D_Ydata = ...
            OD_Smooth(CloudCntr(2),:);
    
        OD_1D_Yfit = ...
            roiDistImage(CloudCntr(2),:);
    
    end


    %% ============================================================
    % 8G. PLOT 1D CROSS-SECTIONS
    % =============================================================

    set( ...
        0, ...
        'CurrentFigure', ...
        figStruct.fig1DFit);


    fit1DAxH(j) = ...
        subplot(nRows,nCols,j);


    hold on;


    % Experimental data.

    plot( ...
        OD_1D_Xdata, ...
        'k.','MarkerSize', 15);

    plot( ...
        OD_1D_Ydata, ...
        'r.','MarkerSize', 15);


    % Fitted distributions.

    plot( ...
        OD_1D_Xfit, ...
        'k');

    plot( ...
        OD_1D_Yfit, ...
        'r');


    if isfield(analyVar,'dimenReduc1D') && ...
            analyVar.dimenReduc1D == 1
    
        title( ...
            sprintf( ...
                '%s, integrated X/Y profiles', ...
                plotLabel), ...
            'Interpreter','none');
            xlabel('Pixel');
            ylabel('Integrated OD');
    
    else
    
        title( ...
            sprintf( ...
                '%s, center [%d,%d]', ...
                plotLabel, ...
                CloudCntr(1), ...
                CloudCntr(2)), ...
            'Interpreter','none');
        xlabel('Pixel');
        ylabel('OD');
    
    end

    grid off;
    box on;


    %% ------------------------------------------------------------
    % Automatically determine useful Y limits
    % -------------------------------------------------------------

    yValues = [
        OD_1D_Xdata(:)
        OD_1D_Ydata(:)
        OD_1D_Xfit(:)
        OD_1D_Yfit(:)
        ];


    validY = ...
        yValues(isfinite(yValues));


    if ~isempty(validY)

        yPadding = ...
            0.05 * max( ...
                max(validY)-min(validY), ...
                1e-6);


        ylim([
            min(validY)-yPadding
            max(validY)+yPadding
            ]);

    end


    hold off;

end


%% ================================================================
%  9. FIGURE TITLES
%  ================================================================

%% Averaged OD image figure

set( ...
    0, ...
    'CurrentFigure', ...
    figStruct.avgODImages);


set( ...
    gcf, ...
    'Name', ...
    'Averaged OD Images: Scan Parameters');


mtit( ...
    'Averaged OD Images: Scan Parameters', ...
    'FontSize',16, ...
    'zoff',.025, ...
    'xoff',-.01);


%% 2D fit figure

set( ...
    0, ...
    'CurrentFigure', ...
    figStruct.fig2DFit);


set( ...
    gcf, ...
    'Name', ...
    '2D Cloud Fit: Averaged Scans');


mtit( ...
    '2D Cloud Fit: Averaged Scans', ...
    'FontSize',16, ...
    'zoff',.025, ...
    'xoff',-.01);


%% Residual figure

set( ...
    0, ...
    'CurrentFigure', ...
    figStruct.figRes);


set( ...
    gcf, ...
    'Name', ...
    'Residuals: Averaged Scans');


mtit( ...
    'Residuals: Averaged Scans', ...
    'FontSize',16, ...
    'zoff',.025, ...
    'xoff',-.01);


%% 1D fit figure

set( ...
    0, ...
    'CurrentFigure', ...
    figStruct.fig1DFit);


oneDTitle = [
    'Cross-Section of Averaged Fit using ' ...
    strrep( ...
        analyVar.fitModel, ...
        analyVar.InitCase, ...
        [analyVar.InitCase ' '])
    ];


set( ...
    gcf, ...
    'Name', ...
    '1D Fit: Averaged Scans');


mtit( ...
    oneDTitle, ...
    'FontSize',16, ...
    'zoff',.025, ...
    'xoff',-.01);


end



%% =================================================================
%  GET DATA USED FOR PLOTTING
%  =================================================================

function plotData = ...
    get_averaged_fit_plot_data( ...
        analyVar,avgDataset)
% Collect averaged OD images, fit parameters, labels, and provenance.


plotData = struct;


%% Make sure the required averaged data exists

if ~isfield(avgDataset,'All_OD_Image')

    error( ...
        'avgDataset.All_OD_Image is missing.');

end


if ~isfield(avgDataset,'All_PCell')

    error( ...
        'avgDataset.All_PCell is missing.');

end


%% Copy the image and fit data

plotData.images = ...
    avgDataset.All_OD_Image;

plotData.PCell = ...
    avgDataset.All_PCell;


%% Determine number of plots

plotData.nPlots = min( ...
    numel(plotData.images), ...
    numel(plotData.PCell));


if plotData.nPlots == 0

    error( ...
        'avgDataset contains no averaged images or fit parameters.');

end


%% Allocate provenance information

plotData.referenceBasenameNums = ...
    ones(plotData.nPlots,1);

plotData.referenceImageNums = ...
    ones(plotData.nPlots,1);

plotData.labels = ...
    cell(plotData.nPlots,1);


%% Build labels for each averaged scan point

for j = 1:plotData.nPlots


    %% Recover source information when available

    if isfield(avgDataset,'sourceInfo') && ...
            numel(avgDataset.sourceInfo) >= j && ...
            ~isempty(avgDataset.sourceInfo{j})

        src = ...
            avgDataset.sourceInfo{j};


        if isfield(src,'basenameNums') && ...
                ~isempty(src.basenameNums)

            plotData.referenceBasenameNums(j) = ...
                src.basenameNums(1);

        end


        if isfield(src,'imageNums') && ...
                ~isempty(src.imageNums)

            plotData.referenceImageNums(j) = ...
                src.imageNums(1);

        end

    end


    %% Determine scan parameter value

    if isfield(avgDataset,'paramVals') && ...
            numel(avgDataset.paramVals) >= j

        parameterValue = ...
            avgDataset.paramVals(j);

    elseif isfield(avgDataset,'imagevcoAtom') && ...
            numel(avgDataset.imagevcoAtom) >= j

        parameterValue = ...
            avgDataset.imagevcoAtom(j);

    else

        parameterValue = j;

    end


    %% Determine tweezer number for display

    if isfield(avgDataset,'sourceInfo') && ...
            numel(avgDataset.sourceInfo) >= j && ...
            ~isempty(avgDataset.sourceInfo{j}) && ...
            isfield( ...
                avgDataset.sourceInfo{j}, ...
                'tweezerNums') && ...
            ~isempty( ...
                avgDataset.sourceInfo{j}.tweezerNums)

        tweezerNum = ...
            avgDataset.sourceInfo{j}.tweezerNums(1);

    elseif isfield(avgDataset,'numTweezers') && ...
            avgDataset.numTweezers > 0

        tweezerNum = ...
            mod(j-1,avgDataset.numTweezers)+1;

    else

        tweezerNum = 1;

    end


    %% Create human-readable plot label

    plotData.labels{j} = ...
        sprintf( ...
            'Tweezer %d, %s = %g %s', ...
            tweezerNum, ...
            analyVar.avgScanParam, ...
            parameterValue, ...
            analyVar.xDataUnit);

end


%% Storage for reconstructed fit windows

plotData.fitWindowCells = ...
    cell(plotData.nPlots,1);

end



%% =================================================================
%  RECONSTRUCT 2D FIT IMAGE
%  =================================================================

function [roiDistImage,OD_Fit_ImageCell] = ...
    reconstruct_fit_image( ...
        analyVar,OD_Image,PCell,basenameNum)
% Re-evaluate the saved fit parameters on the image grid.


OD_Image = ...
    double(OD_Image);


%% Apply the same smoothing used during fitting

if analyVar.fitSmoothOD

    OD_FitInput = ...
        analyVar.smoothFilt( ...
            OD_Image, ...
            analyVar.smoothFiltMat);

else

    OD_FitInput = ...
        OD_Image;

end


%% Tweezer image is already the complete fit window

if isfield(analyVar,'UseTweezer') && ...
        analyVar.UseTweezer == 1

    fitWinSize = ...
        size(OD_Image,1);


    if size(OD_Image,2) ~= fitWinSize

        error( ...
            'Tweezer OD image must be square.');

    end


    OD_Fit_ImageCell = {
        OD_FitInput
        };


    roiWin_Index = {
        true(size(OD_Image))
        };

else

    error([ ...
        'The current averaged-fit plotting helper expects ' ...
        'tweezer images with UseTweezer = 1.']);

end


%% Calculate weighting/error image

errCell = ...
    cellfun( ...
        @(fitParameters,fitImage) ...
            get_OD_weight( ...
                fitParameters(end), ...
                fitImage), ...
        PCell, ...
        OD_Fit_ImageCell, ...
        'UniformOutput',false);


%% Build model coordinate grid

[Xgrid,Ygrid] = ...
    meshgrid(1:fitWinSize);


%% Evaluate the saved fit parameters

fitDistCell = ...
    cellfun( ...
        @(fitParameters,fitError) ...
            reshape( ...
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

    roiDistImage = ...
        fitDistCell{1};

else

    roiDistImage = ...
        zeros(size(OD_Image));


    for fitNum = 1:numel(fitDistCell)

        roiMask = ...
            roiWin_Index{fitNum};

        roiMask(roiMask) = ...
            fitDistCell{fitNum};

        roiDistImage = ...
            roiDistImage + roiMask;

    end

end

end



%% =================================================================
%  GET FITTED CLOUD CENTER
%  =================================================================

function CloudCntr = ...
    get_fit_center_from_pcell( ...
        analyVar, ...
        PCell, ...
        imageSize, ...
        basenameNum, ...
        isLocalTweezerImage)
% Extract x/y center parameters from the fitted parameter cell.


fitParameters = ...
    PCell{1};


parameterNames = ...
    analyVar.InitCondList;


if isstring(parameterNames)

    parameterNames = ...
        cellstr(parameterNames);

end


%% Find X center parameter

xIndex = ...
    find( ...
        strcmpi(parameterNames,'xCntr') | ...
        strcmpi(parameterNames,'xCenter') | ...
        strcmpi(parameterNames,'x0'), ...
        1);


%% Find Y center parameter

yIndex = ...
    find( ...
        strcmpi(parameterNames,'yCntr') | ...
        strcmpi(parameterNames,'yCenter') | ...
        strcmpi(parameterNames,'y0'), ...
        1);


%% Fall back to image center if parameters are unavailable

if isempty(xIndex) || isempty(yIndex)

    CloudCntr = [
        (imageSize(2)+1)/2
        (imageSize(1)+1)/2
        ];

    return;

end


%% Extract fitted center

fitX = ...
    fitParameters(xIndex);

fitY = ...
    fitParameters(yIndex);


%% Convert coordinates if necessary

if isLocalTweezerImage

    % The fit coordinates already correspond to the local
    % tweezer image.

    CloudCntr = [
        fitX
        fitY
        ];

else

    % Preserve the original full-ROI coordinate convention.

    CloudCntr = ...
        ( ...
        analyVar.roiWinRadAtom(basenameNum) - ...
        (analyVar.funcFitWin(basenameNum)-1)/2 ...
        ) + ...
        [fitX fitY];

end

end



%% =================================================================
%  GLOBAL IMAGE COLOR LIMITS
%  =================================================================

function colorLimits = ...
    get_global_image_limits(imageCell)
% Find the minimum and maximum finite values across all images.


globalMin = inf;
globalMax = -inf;


for j = 1:numel(imageCell)

    imageData = ...
        imageCell{j};


    if isempty(imageData)

        continue;

    end


    validData = ...
        imageData(isfinite(imageData));


    if isempty(validData)

        continue;

    end


    globalMin = ...
        min(globalMin,min(validData));

    globalMax = ...
        max(globalMax,max(validData));

end


%% Return valid limits

if isfinite(globalMin) && ...
        isfinite(globalMax) && ...
        globalMin < globalMax

    colorLimits = [
        globalMin
        globalMax
        ];

else

    colorLimits = [
        NaN
        NaN
        ];

end

end



%% =================================================================
%  SYMMETRIC RESIDUAL COLOR LIMITS
%  =================================================================

function colorLimits = ...
    get_symmetric_global_image_limits(imageCell)
% Find the largest absolute residual and use +/- that value.


largestMagnitude = 0;


for j = 1:numel(imageCell)

    imageData = ...
        imageCell{j};


    if isempty(imageData)

        continue;

    end


    validData = ...
        imageData(isfinite(imageData));


    if isempty(validData)

        continue;

    end


    largestMagnitude = ...
        max( ...
            largestMagnitude, ...
            max(abs(validData)));

end


%% Return symmetric limits

if largestMagnitude > 0

    colorLimits = [
        -largestMagnitude
         largestMagnitude
        ];

else

    colorLimits = [
        NaN
        NaN
        ];

end

end