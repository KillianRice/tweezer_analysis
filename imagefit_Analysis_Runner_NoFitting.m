% Script to run complete analysis of Neutral table data
% 
% Please be sure to have the correct date defined in AnalysisVariables
% prior to running.

% clear current figures
close all

% Load variables and file data
analyVar = AnalysisVariables;
if analyVar.numBasenamesAtom > 5
    disp(['About to analyze ' num2str(analyVar.numBasenamesAtom) ' scans, are you sure you want to continue?']);
    s = input('Y/N >','s');
    while s ~= 'Y' && s ~= 'N' && s ~= 'y' && s ~= 'n'
        s = input('Y/N >','s');
        disp(s)
    end
    if s == 'n' || s == 'N'
        return;
    end
end
indivDataset = get_indiv_batch_data(analyVar);

AvgOfroiImageAtom = repmat({0},1,analyVar.numBasenamesAtom);
AvgOfBackground = repmat({0},1,analyVar.numBasenamesAtom);
AvgOfroiImageAtomsMinusBkg = repmat({0},1,analyVar.numBasenamesAtom);
AvgOfroiImageAtomsMinusBkgErr = repmat({0},1,analyVar.numBasenamesAtom);


SumAvgOfroiImageAtom = repmat({0},1,analyVar.numBasenamesAtom);
SumAvgOfBackground = repmat({0},1,analyVar.numBasenamesAtom);
SumAvgOfroiImageAtomsMinusBkg = repmat({0},1,analyVar.numBasenamesAtom);

for basenameNum = 1:analyVar.numBasenamesAtom % number of files to analyze.
    SumOfroiImageAtom = 0;
    SumOfBackground= 0;
    SumOfroiImageAtomsMinusBkg = 0;
    SumSqOfroiImageAtomsMinusBkg = 0;
    for k = 1:indivDataset{basenameNum}.CounterAtom % number of experimental shots in each file

        % Grabs an ROI cut of each of the images
        roiImageAtom = reshape(indivDataset{basenameNum}.AtomsCloud(:,k),[1,1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1)); 
        roiImageBackground = reshape(indivDataset{basenameNum}.BackgroundCloud(:,k),[1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1)); 
        
        % Add each image into a summed total
        SumOfroiImageAtom = SumOfroiImageAtom + roiImageAtom;
        SumOfroiImageAtomsMinusBkg = SumOfroiImageAtomsMinusBkg + (roiImageAtom - roiImageBackground);
        SumOfBackground = SumOfBackground + roiImageBackground;
        SumSqOfroiImageAtomsMinusBkg = SumSqOfroiImageAtomsMinusBkg + (roiImageAtom - roiImageBackground).^2;

        %%% Save OD
        %dlmwrite([analyVar.analyOutDir char(indivDataset{basenameNum}.fileAtom(k)) analyVar.ODimageFilename],(roiImageAtom - roiImageBackground),'\t');

    end
    ADU2electron = 0.45;

    N = indivDataset{basenameNum}.CounterAtom;
    
    AvgOfroiImageAtom{1,basenameNum} = ADU2electron * SumOfroiImageAtom / N;
    AvgOfBackground{1,basenameNum} = ADU2electron * SumOfBackground / N;
    AvgOfroiImageAtomsMinusBkg{1,basenameNum} = ADU2electron * SumOfroiImageAtomsMinusBkg / N;
    
    SumAvgOfroiImageAtom{1,basenameNum} = sum(AvgOfroiImageAtom{1,basenameNum}(:));
    SumAvgOfBackground{1,basenameNum} = sum(AvgOfBackground{1,basenameNum}(:));
    SumAvgOfroiImageAtomsMinusBkg{1,basenameNum} = sum(AvgOfroiImageAtomsMinusBkg{1,basenameNum}(:));

    %finding error for each pixel of each array
    
    % Variance per pixel
    N = indivDataset{basenameNum}.CounterAtom;
    varImage = (ADU2electron*SumSqOfroiImageAtomsMinusBkg - ((ADU2electron*SumOfroiImageAtomsMinusBkg).^2)/N) / (N - 1);
    
    % Standard error of the mean per pixel
    semImage = sqrt(varImage / N);
    AvgOfroiImageAtomsMinusBkgErr{1,basenameNum} = sqrt(sum(semImage(:).^2));
end


dataCells = {AvgOfroiImageAtom, AvgOfBackground, AvgOfroiImageAtomsMinusBkg};  % cell array of 2D matrices
dataNames = {'Atom', 'Background', 'Atom - Background'};
numPerFig = 10;              % 5x2 layout
n = numel(dataCells{1});


for i = 1:n
    % Start a new figure every 15 plots
    if mod(i-1, numPerFig) == 0
        figure;
        t = tiledlayout(5,2, 'Padding', 'compact', 'TileSpacing', 'compact');
        name = "Histograms of File" + i + " to " + min(i+9,n);
        title(t, sprintf(name));
    end

    % Get current matrix
    A = dataCells{1}{i};
    B = dataCells{2}{i};
    C = dataCells{3}{i};

    %% Plot Individual histogram

    % Select next tile
    %nexttile;
    %histogram(A(:), 100, 'BinLimits', [0, 250]);
    %title(sprintf('Atom Array %d', i));
    %% Select next tile
    %nexttile;
    %histogram(B(:), 100, 'BinLimits', [0, 250]);
    %title(sprintf('Background Array %d', i));
    %% Select next tile
    %nexttile;
    %histogram(C(:), 100, 'BinLimits', [-50, 50]);

    %% Keep All 3 on same histogram
    nexttile;
    hold on;
    edges = linspace(-5, 100, 400);
    histogram(A, edges, 'FaceColor', 'r', 'FaceAlpha', 0.3);
    histogram(B, edges, 'FaceColor', 'b', 'FaceAlpha', 0.3);
    histogram(C, edges, 'FaceColor', 'g', 'FaceAlpha', 0.3);
    hold off;
    xlabel('e- Counts')
    ylabel('Number of Pixels')
    legend('Atom Img', 'Bkg Img', 'Atom - Bkg Img');
    yscale log
    title(sprintf('File %d', analyVar.meanListVar(i)));


end

%% Plot all scans on top of each other for each category
plotAll = 1;
if plotAll
    figure;
    t = tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');
    name = "Histograms of Files: Atom, Bkg, Atom - Bkg";
    title(t, sprintf(name));
    
    nexttile;
    hold on;
    edges = linspace(40, 52, 100);
    for i = 1:n
        histogram(dataCells{1}{i}, edges, 'FaceColor', analyVar.COLORS(i,:), 'FaceAlpha', 0.3);
    end
    hold off;
    xlabel('e- Counts')
    ylabel('Number of Pixels')
    yscale log
    legend('532 On', '532 Off');
    title(sprintf('Atoms'));

    nexttile;
    hold on;
    edges = linspace(40, 50, 100);
    for i = 1:n
        histogram(dataCells{2}{i}, edges, 'FaceColor', analyVar.COLORS(i,:), 'FaceAlpha', 0.3);
    end
    hold off;
    xlabel('e- Counts')
    ylabel('Number of Pixels')
    yscale log
    legend('532 On', '532 Off');
    title(sprintf('Bkg'));

    nexttile;
    hold on;
    edges = linspace(-1, 6, 100);
    for i = 1:n
        histogram(dataCells{3}{i}, edges, 'FaceColor', analyVar.COLORS(i,:), 'FaceAlpha', 0.3);
    end
    hold off;
    xlabel('e- Counts')
    ylabel('Number of Pixels')
    yscale log
    legend('532 On', '532 Off');
    title(sprintf('Atoms - Bkg'));


end





%Subtract the atoms minus the background for plotting and prepare an array
%to plot transects of the image
numPerFig = 5;
for i = 1:n

    
    % Start a new figure every 15 plots
    if mod(i-1, numPerFig) == 0
            figure;
            t = tiledlayout(5,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    end
    nexttile
    imshow(AvgOfroiImageAtom{1,i}, [], 'Colormap', jet)
    h = gca;
    h.Visible = 'On';
    colorbar
    title("Atom Image")

    nexttile
    imshow(AvgOfBackground{1,i}, [], 'Colormap', jet)
    colorbar
    title("Background Image")
    h = gca;
    h.Visible = 'On';

    nexttile
    imshow(AvgOfroiImageAtomsMinusBkg{1,i}, [], 'Colormap', hot)
    colorbar
    h = gca;
    h.Visible = 'On';
    title("Atom - Background Image")
end


scanIDs = analyVar.meanListVar;
[sorted_ID, idx] = sort(scanIDs);

%% Atom Image
figure;
listOfTotalADU = cell2mat(SumAvgOfroiImageAtom);
sorted_TotalADU = listOfTotalADU(idx);
listOfTotalADUErr = cell2mat(AvgOfroiImageAtomsMinusBkgErr);
sorted_TotalADUErr = listOfTotalADUErr(idx);
errorbar(sorted_ID, sorted_TotalADU, sorted_TotalADUErr,...
        'LineStyle','-',...
        'Marker', analyVar.MARKERS2(1),...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(1,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(1,:));
title('Avg Atom Img e- on Zyla');
xlabel('RMOT X axis');
ylabel('Avg Total Number of e-');
hold off

%% Background Image
figure;
listOfTotalADU = cell2mat(SumAvgOfBackground);
sorted_TotalADU = listOfTotalADU(idx);
listOfTotalADUErr = cell2mat(AvgOfroiImageAtomsMinusBkgErr);
sorted_TotalADUErr = listOfTotalADUErr(idx);
errorbar(sorted_ID, sorted_TotalADU, sorted_TotalADUErr,...
        'LineStyle','-',...
        'Marker', analyVar.MARKERS2(1),...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(1,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(1,:));
title('Avg Bkg e- on Zyla');
xlabel('Avg Rep');
ylabel('Avg Total Number of e-');
hold off


%% Atom - Background Image
figure;
listOfTotalADU = cell2mat(SumAvgOfroiImageAtomsMinusBkg);
sorted_TotalADU = listOfTotalADU(idx);
listOfTotalADUErr = cell2mat(AvgOfroiImageAtomsMinusBkgErr);
sorted_TotalADUErr = listOfTotalADUErr(idx);
errorbar(sorted_ID, sorted_TotalADU, sorted_TotalADUErr,...
        'LineStyle','-',...
        'Marker', analyVar.MARKERS2(1),...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(1,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(1,:));
title('Avg Atom - Bkg e- on Zyla');
xlabel('Avg Rep');
ylabel('Avg Total Number of e-');
hold off


%plot a transect based on where the user clicks on the atoms-background
%image

% sizeOfImage = size(atomsMinusBackground);
% xTransectArray= linspace(0, sizeOfImage(1) - 1, sizeOfImage(1));
% yTransectArray= linspace(0, sizeOfImage(2) - 1, sizeOfImage(2));
% 
% [xTransect,yTransect] = getpts(gca)
% 
% xTransect = round(xTransect);
% yTransect = round(yTransect);
% 
% xPntsForLine = zeros(sizeOfImage(1),1) + xTransect;
% yPntsForLine = zeros(sizeOfImage(2),1) + yTransect;
% 
% plot(xTransectArray, yPntsForLine)
% plot(xPntsForLine, yTransectArray)
% 
% 
% nexttile
% plot(xTransectArray,atomsMinusBackground(yTransect, :))
% xlabel = ('x position')
% ylabel = ('Proportional to Counts')
% title("X-transect of Atom - Background at ", xTransect)
% 
% nexttile
% plot(yTransectArray,atomsMinusBackground(:, xTransect))
% xlabel = ('y position')
% ylabel = ('Proportional to Counts')
% title("Y-transect of Atom - Background at ", yTransect)


%imshow(SumOfroiImageNotAtom,[])
% % Plotting routine
% if analyVar.SavePlotData == 1
%     % Save data from output
%     PlotData = imagefit_ParamEval(analyVar,indivDataset);
% else
%     imagefit_ParamEval(analyVar,indivDataset);
%end

%fitBackgroundGaussian(dataCells{1},400)
%fitBackgroundGaussian(dataCells{2},400)
%fitBackgroundGaussian(dataCells{3},400)

function params = fitBackgroundGaussian(imageList, nBins)

if nargin < 2
    nBins = 100;
end

nImages = length(imageList);

% Preallocate struct array
params(nImages) = struct('A', [], 'mu', [], 'sigma', []);

for i = 1:nImages
    
    img = imageList{i};
    pixels = img(:);
    pixels = pixels(isfinite(pixels));
    
    % Histogram
    [counts, edges] = histcounts(pixels, nBins);
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;

    % Gaussian model
    gaussEqn = fittype('A*exp(-(x-mu)^2/(2*sigma^2))', ...
        'independent', 'x', ...
        'coefficients', {'A', 'mu', 'sigma'});

    % Initial guesses
    A0 = max(counts);
    mu0 = mean(pixels);
    sigma0 = std(pixels);

    % Fit
    fitResult = fit(binCenters(:), counts(:), gaussEqn, ...
        'StartPoint', [A0, mu0, sigma0]);

    % Store parameters
    params(i).A = fitResult.A;
    params(i).mu = fitResult.mu;
    params(i).sigma = fitResult.sigma;

    % Plot (one figure per image)
    figure;
    bar(binCenters, counts, 'hist');
    hold on;

    xFit = linspace(min(binCenters), max(binCenters), 500);
    yFit = fitResult.A * exp(-(xFit - fitResult.mu).^2 / (2 * fitResult.sigma^2));

    plot(xFit, yFit, 'r-', 'LineWidth', 2);

    xlabel('Signal (e-)');
    ylabel('Counts');
    title(sprintf('Image %d: Histogram + Gaussian Fit', i));

    legend('Data', sprintf('\\mu = %.3f, \\sigma = %.3f', ...
        fitResult.mu, fitResult.sigma));

    grid on;
    hold off;

end

end

results = analyzeHistCell_EM(dataCells{3}, true);
plot(results.LLR, '-o');
title('Detection significance vs parameter');
function results = analyzeHistCell_EM(cellData, plotFlag)

nCells = numel(cellData);

results.lambda0 = zeros(nCells,1);
results.lambda1 = zeros(nCells,1);
results.weight  = zeros(nCells,1);
results.snr     = zeros(nCells,1);
results.LLR     = zeros(nCells,1); % likelihood ratio

for k = 1:nCells
    
    data = cellData{k}(:);
    % for j = 1:numel(data)
    %     if data(j) < 0
    %         data(j) = 0;
    %     end
    % end
    data = data + abs(min(data));
    data = data(~isnan(data));
    
    % --- INITIALIZATION ---
    lambda0 = mean(data)*0.9;
    lambda1 = mean(data)*1.1;
    w = 0.5;
    
    % --- EM ALGORITHM ---
    for iter = 1:200
        
        % E-step
        P0 = poisspdf(data, lambda0);
        P1 = poisspdf(data, lambda1);
        
        gamma = w*P1 ./ (w*P1 + (1-w)*P0 + 1e-12);
        
        % M-step
        w = mean(gamma);
        lambda1 = sum(gamma .* data) / sum(gamma);
        lambda0 = sum((1-gamma).*data) / sum(1-gamma);
        
        % Enforce separation
        if lambda0 > lambda1
            [lambda0, lambda1] = deal(lambda1, lambda0);
            w = 1 - w;
        end
    end
    
    % --- STORE ---
    results.lambda0(k) = lambda0;
    results.lambda1(k) = lambda1;
    results.weight(k)  = w;
    
    % --- SNR (still useful but secondary now) ---
    results.snr(k) = (lambda1 - lambda0) / sqrt(lambda0 + lambda1);
    
    % --- LIKELIHOODS ---
    logL_mix = sum(log(w*poisspdf(data,lambda1) + ...
                       (1-w)*poisspdf(data,lambda0) + 1e-12));
    
    % Single Poisson fit
    lambda_single = mean(data);
    logL_single = sum(log(poisspdf(data, lambda_single) + 1e-12));
    
    % Likelihood Ratio (KEY METRIC)
    results.LLR(k) = 2*(logL_mix - logL_single);
    
    % --- PLOTTING ---
    if plotFlag
        
        disp("Lambdas")
        disp(lambda0)
        disp(lambda1)
        figure; hold on;
        
        edges = linspace(0, 20, 100);
        histogram(data, edges, 'Normalization','pdf', 'FaceAlpha',0.4);
        
        nVals = 0:max(data);
        fit_mix = w*poisspdf(nVals,lambda1) + ...
                  (1-w)*poisspdf(nVals,lambda0);
        
        fit_single = poisspdf(nVals, lambda_single);
        
        plot(nVals, fit_mix, 'r-', 'LineWidth',2);
        plot(nVals, fit_single, 'k--', 'LineWidth',1.5);
        
        title(sprintf(['Cell %d\nSNR=%.2f | LLR=%.2f'], ...
              k, results.snr(k), results.LLR(k)));
        
        legend('Data','Mixture','Single Poisson');
        hold off;
    end
end
end