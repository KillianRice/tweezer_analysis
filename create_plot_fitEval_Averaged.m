function [fit2DAxH, fit1DAxH, resAxH] = create_plot_fitEval_Averaged(analyVar, avgDataset)

%% Figure numbers
figStruct.fig2DFit = analyVar.figNum.fig2DFit;
figStruct.figRes   = analyVar.figNum.figRes;
figStruct.fig1DFit = analyVar.figNum.fig1DFit;
figStruct.fig1DBEC = analyVar.figNum.fig1DBEC;

figure(figStruct.fig2DFit);
figure(figStruct.figRes);
figure(figStruct.fig1DFit);

if strcmpi(analyVar.InitCase,'Bimodal')
    figure(figStruct.fig1DBEC);
end

%% Plot layout
nPlots = avgDataset.CounterAtom;

if isfield(avgDataset,'SubPlotRows') && isfield(avgDataset,'SubPlotCols')
    nRows = avgDataset.SubPlotRows;
    nCols = avgDataset.SubPlotCols;
else
    nRows = ceil(sqrt(nPlots));
    nCols = ceil(nPlots/nRows);
end

%% Preallocate axes
[fit2DAxH, fit1DAxH, resAxH] = deal(zeros(1,nPlots));

%% Choose reference basename index for window sizes
basenameNum = 1;

for j = 1:nPlots

    %% Preallocate reconstructed fit image
    roiDistImage = zeros(size(avgDataset.All_OD_Image{j}));

    %% Smooth averaged OD image
    OD_Smooth = analyVar.smoothFilt(avgDataset.All_OD_Image{j}, analyVar.smoothFiltMat);

    %% Get ROI windows
    if isfield(avgDataset,'roiWin_Index')
        roiWin_Index = avgDataset.roiWin_Index;
    else
        error('avgDataset.roiWin_Index is missing. Save roiWin_Index into avgDataset during averaging.');
    end

    %% Retrieve fitted OD windows
    if isfield(analyVar,'UseTweezer') && analyVar.UseTweezer == 1

        OD_Image_Single = dlmread(avgDataset.avgODFiles{j});
    
        roiWin_Index = {true(size(OD_Image_Single))};
        fitWinSize = size(OD_Image_Single,1);

        OD_Fit_ImageCell = cellfun( ...
            @(x) reshape(OD_Image_Single(x), [1 1]*fitWinSize), ...
            roiWin_Index, ...
            'UniformOutput', 0);
    
    else
        OD_Fit_ImageCell = cellfun( ...
            @(x) reshape(avgDataset.All_OD_Image{j}(x), [1 1]*analyVar.funcFitWin(basenameNum)), ...
            cellfun(analyVar.fitWinLogicInd, roiWin_Index, 'UniformOutput', 0), ...
            'UniformOutput', 0);
    end

    %% Retrieve fit coefficients
    PCell = avgDataset.All_PCell{j};

    %% Weighting
    errCell = cellfun(@(x,y) get_OD_weight(x(end), y), ...
        PCell, OD_Fit_ImageCell, 'UniformOutput', 0);

    %% Generate fit model image
    [Xgrid,Ygrid] = meshgrid(1:fitWinSize);%analyVar.funcFitWin(basenameNum));

    fitDistCell = cellfun( ...
        @(x,y) reshape( ...
            feval(str2func(analyVar.fitModel), ...
            x(1:length(analyVar.InitCondList)), ...
            [Xgrid(:), Ygrid(:), y(:)]), ...
            [1 1]*fitWinSize), ... %analyVar.funcFitWin(basenameNum
        PCell, errCell, 'UniformOutput', 0);

    %% Insert each fit window into full ROI image
    for i = 1:sum(analyVar.LatticeAxesFit)

        roiDistTmp = roiWin_Index{i};

        roiDistTmp(analyVar.fitWinLogicInd(roiDistTmp)) = fitDistCell{i};

        roiDistImage = roiDistImage + roiDistTmp;
    end

    %% Label for this averaged point
    if isfield(avgDataset,'paramVals')
        plotLabel = sprintf('%g %s', avgDataset.paramVals(j), analyVar.xDataUnit);
    else
        plotLabel = sprintf('%s %g',analyVar.avgScanParam, avgDataset.imagevcoAtom(j));
    end

    %% Plot 2D fit
    set(0,'CurrentFigure',figStruct.fig2DFit)
    fit2DAxH(j) = subplot(nRows,nCols,j);

    pcolor(roiDistImage);
    shading flat;
    colorbar;
    title(plotLabel);
    grid off;

    if j == nPlots
        set(gcf,'Name','2D Cloud Fit: Averaged Scans');
        mtit('2D Cloud Fit: Averaged Scans','FontSize',16,'zoff',.025,'xoff',-.01);
    end

    %% Plot residuals
    set(0,'CurrentFigure',figStruct.figRes)
    resAxH(j) = subplot(nRows,nCols,j);

    pcolor(avgDataset.All_OD_Image{j} - roiDistImage);
    shading flat;
    colorbar;
    title(plotLabel);
    grid off;

    if j == nPlots
        set(gcf,'Name','Residuals: Averaged Scans');
        mtit('Residuals: Averaged Scans','FontSize',16,'zoff',.025,'xoff',-.01);
    end

    %% Plot 1D cross sections
    set(0,'CurrentFigure',figStruct.fig1DFit)

    CloudCntr = round( ...
        (analyVar.roiWinRadAtom(basenameNum) - ...
        (analyVar.funcFitWin(basenameNum) - 1)/2) + ...
        [avgDataset.All_fitParams{j}{1}.xCntr, ...
         avgDataset.All_fitParams{j}{1}.yCntr]);

    % Clamp center to image bounds
    CloudCntr(1) = max(1, min(size(avgDataset.All_OD_Image{j},2), CloudCntr(1)));
    CloudCntr(2) = max(1, min(size(avgDataset.All_OD_Image{j},1), CloudCntr(2)));

    OD_1D_Xdata = OD_Smooth(:,CloudCntr(1));
    OD_1D_Ydata = OD_Smooth(CloudCntr(2),:);

    OD_1D_Xfit = roiDistImage(:,CloudCntr(1));
    OD_1D_Yfit = roiDistImage(CloudCntr(2),:);

    fit1DAxH(j) = subplot(nRows,nCols,j);

    hold on;
    grid off;

    plot(OD_1D_Xdata,'c.');
    plot(OD_1D_Ydata,'g.');
    plot(OD_1D_Xfit,'k');
    plot(OD_1D_Yfit,'r');

    title(sprintf('%s [%g,%g]', plotLabel, CloudCntr(1), CloudCntr(2)));

    xlim(mean(CloudCntr) + [-1 1]*analyVar.roiWinRadAtom(basenameNum));

    ylim([-0.1, ...
        max(max([OD_1D_Xdata; OD_1D_Ydata'; OD_1D_Xfit; OD_1D_Yfit'])) + 0.1]);

    if j == nPlots
        legend('y data','x data','y fit','x fit');
        set(gcf,'Name','1D Fit: Averaged Scans');

        mtit(['Cross-Section of Averaged Fit using ' ...
            strrep(analyVar.fitModel, analyVar.InitCase, [analyVar.InitCase ' '])], ...
            'FontSize',16,'zoff',.025,'xoff',-.01);
    end

end

end