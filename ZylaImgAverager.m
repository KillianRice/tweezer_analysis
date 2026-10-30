
%% Use for grabbing the raw csv file that holds all the images to average into single image
% % 1. Read the CSV file into a matrix
% % For very large files, readmatrix is highly optimized
% fprintf('Reading CSV file... \n');
% data = readmatrix('ZylaUnblocked461OnImages.csv');
% 
% % 2. Average all columns into a single column (dimension 2)
% fprintf('Averaging columns... \n');
% columnAverage = mean(data, 2); 
% 
% % Clear the massive raw data matrix to free up system memory
% clear data; 
% 
% % 3. Calculate 'n' dynamically for the nxn grid
% totalRows = length(columnAverage);
% n = sqrt(totalRows);
% 
% if rem(n, 1) ~= 0
%     error('The number of rows (%d) cannot form a perfect square (nxn) image.', totalRows);
% end
% 
% % 4. Reshape the mx1 array into an nxn matrix
% % Note: MATLAB reshapes column-wise by default. 
% % We use .' (transpose) if you want the rows filled sequentially from left-to-right.
% fprintf('Reshaping into a %dx%d image... \n', n, n);
% imageMatrix = reshape(columnAverage, n, n).'; 
% 
% % 5. Plot the array as an image
% figure;
% imagesc(imageMatrix);            % Display the scaled data as an image
% colormap('gray');                % Set color scale (e.g., 'gray', 'jet', 'viridis')
% colorbar;                        % Show the color intensity scale bar
% axis image;                      % Ensure pixel aspect ratio is 1:1
% title(sprintf('Visualized %dx%d Data Grid', n, n));
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');

%% Use for grabbing the average image to calculate the std for each pixel
% % 1. Read the combined CSV file
% fprintf('Reading CSV file... \n');
% data = readmatrix('ZylaUnblocked461OnImages.csv');
% 
% % 2. Calculate the Standard Deviation (STD) for each row (dimension 2)
% fprintf('Calculating row-wise standard deviations... \n');
% rowSTD = std(data, 0, 2); % '0' specifies the default normalization (N-1)
% 
% % Clear the massive raw data matrix immediately to free up system memory
% clear data; 
% 
% % 3. Define dimensions (4194304 rows = 2048 x 2048)
% n = 300;
% 
% % 4. Reshape the 4194304x1 STD array into a 2048x2048 matrix
% % Using .' (transpose) matches the row-by-row layout of the previous code
% fprintf('Reshaping into a %dx%d STD image... \n', n, n);
% imageSTD = reshape(rowSTD, n, n).'; 
% 
% writematrix(imageSTD, 'ZylaUnblocked461OnSD.csv');
% 
% % 5. Plot the Standard Deviation Map
% figure;
% imagesc(imageSTD);            
% colormap('hot');                 % 'hot' or 'jet' colormaps are great for variability/uncertainty
% colorbar;                        
% axis image;                      
% title('Pixel Standard Deviation (STD) Map');
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');

%% Use for plotting the averaged image and the 2D array of each pixel's std
% close all;
% 
% fprintf('Reading CSV file... \n');
% avgImage = readmatrix('ZylaBkgSubAvg.csv');
% avgSD = readmatrix('ZylaBkgSubSD.csv');
% 
% % test = readmatrix('ZylaBkgSubTest600.csv');
% % avgImage = test - avgImage;
% 
% avgImage = rot90(avgImage);
% avgSD = rot90(avgSD);
% 
% figure;
% imagesc(avgImage);            
% colormap('hot');
% colorbar;                        
% axis image;
% clim([-10 10])
% title('Random Raw Img - Avg Bkg Zyla Pixel Map');
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');
% 
% figure;
% imagesc(avgSD);            
% colormap('hot');
% colorbar;                        
% axis image;
% clim([0 20])
% title('Pixel Standard Deviation (STD) Map');
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');
% 
% fprintf('Min Standard Deviation: %.4f\n', min(avgSD(:)));
% fprintf('Max Standard Deviation: %.4f\n', max(avgSD(:)));
% fprintf('Mean Standard Deviation: %.4f\n', mean(avgSD(:)));
% fprintf('Median Standard Deviation: %.4f\n', median(avgSD(:)));
% 
% fprintf('Min Pixel Count: %.4f\n', min(avgImage(:)));
% fprintf('Max Pixel Count: %.4f\n', max(avgImage(:)));
% fprintf('Mean Pixel Count: %.4f\n', mean(avgImage(:)));
% fprintf('Median Pixel Count: %.4f\n', median(avgImage(:)));
% 
% figure;
% histogram(avgImage(:))
% %set(gca, 'YScale', 'log');
% ylim([0.5, Inf]);
% xlim([-20, 20]);
% title('Random Raw Img - Avg Bkg Pixel Histogram');
% xlabel('Counts');
% ylabel('# Pixels');
% 
% 
% figure;
% histogram(avgSD(:))
% set(gca, 'YScale', 'log');
% ylim([0.5, Inf]);
% xlim([0, 40]);
% title('Pixel Standard Deviation Histogram');
% xlabel('Standard Deviation');
% ylabel('# Pixels');
% 
% std(avgImage(:))


%% Subtract average image from whole data set
% % 1. Read the CSV file into a matrix
% % For very large files, readmatrix is highly optimi
% fprintf('Reading CSV file... \n');
% data = readmatrix('ZylaBkgData994Images.csv');
% avgImage = readmatrix('AvgZylaBkg.csv').'; %% need to untranpose the image since it was saved transpose
% 
% %Define dimensions (4194304 rows = 2048 x 2048)
% n = 2048;
% [numRows, numCols] = size(data);
% 
% % Save some raw individual data
% test = data(:,1);
% test = reshape(test, n, n).';
% writematrix(test, 'ZylaBkgSubTest1.csv');
% 
% test = data(:,100);
% test = reshape(test, n, n).';
% writematrix(test, 'ZylaBkgSubTest100.csv');
% 
% test = data(:,600);
% test = reshape(test, n, n).';
% writematrix(test, 'ZylaBkgSubTest600.csv');
% 
% 
% % 2. Subtract Avg Img from all raw images
% for i = 1:numCols
% 
%     data(:,i) = data(:,i) - avgImage(:);
% 
% end
% 
% % 3A. Average all columns into a single column (dimension 2)
% fprintf('Averaging columns... \n');
% columnAverage = mean(data, 2); 
% 
% % 3B. Calculate the Standard Deviation (STD) for each row (dimension 2)
% fprintf('Calculating row-wise standard deviations... \n');
% rowSTD = std(data, 0, 2); % '0' specifies the default normalization (N-1)
% 
% clear data;
% 
%  % 4. Reshape the mx1 array into an nxn matrix
%  % Note: MATLAB reshapes column-wise by default. 
%  % We use .' (transpose) if you want the rows filled sequentially from left-to-right.
%  fprintf('Reshaping into a %dx%d image... \n', n, n);
%  imageMatrix = reshape(columnAverage, n, n).';
% 
% % 4. Reshape the 4194304x1 STD array into a 2048x2048 matrix
% % Using .' (transpose) matches the row-by-row layout of the previous code
% fprintf('Reshaping into a %dx%d STD image... \n', n, n);
% imageSTD = reshape(rowSTD, n, n).'; 
% 
% writematrix(imageMatrix, 'ZylaBkgSubAvg.csv');
% writematrix(imageSTD, 'ZylaBkgSubSD.csv');
% 
% % 5. Plot the array as an image
% figure;
% imagesc(imageMatrix);            % Display the scaled data as an image
% colormap('gray');                % Set color scale (e.g., 'gray', 'jet', 'viridis')
% colorbar;                        % Show the color intensity scale bar
% axis image;                      % Ensure pixel aspect ratio is 1:1
% title(sprintf('Visualized %dx%d Data Grid', n, n));
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');
% 
% % 5. Plot the Standard Deviation Map
% figure;
% imagesc(imageSTD);            
% colormap('hot');                 % 'hot' or 'jet' colormaps are great for variability/uncertainty
% colorbar;                        
% axis image;                      
% title('Pixel Standard Deviation (STD) Map');
% xlabel('Pixels (Width)');
% ylabel('Pixels (Height)');
 





%% Fit and plot the bell curve of the histogram
% Prepare data

close all;

avgImage = readmatrix('AvgZylaBkg.csv');
test = readmatrix('ZylaBkgSubTest600.csv');
avgImage = test - avgImage;

data = double(avgImage(:));
data = data(isfinite(data));

if numel(data) > 100000
    rng(1);
    data = data(randperm(numel(data), 100000));
end

% Fit distributions
pdNormal = fitdist(data, 'Normal');
pdT = fitdist(data, 'tLocationScale');
if all(data > 0)
    pdLognormal = fitdist(data, 'Lognormal');
end

% Skew-normal fit: parameters [location, scale, shape]
skewpdf = @(x,p) 2./p(2).*normpdf((x-p(1))./p(2)) ...
    .*normcdf(p(3).*(x-p(1))./p(2));
nll = @(p) -sum(log(max(skewpdf(data,p),realmin)));
q = fminsearch(@(q) nll([q(1),exp(q(2)),q(3)]), ...
    [mean(data),log(std(data)),0]);
p = [q(1),exp(q(2)),q(3)];

% Skew-t fit
q0 = [median(data),log(std(data)),0,log(5)];
opts = optimset('Display','off','MaxIter',3000, ...
    'MaxFunEvals',6000,'TolX',1e-7,'TolFun',1e-7);
qhat = fminsearch(@(q) skewt_nll(data,q),q0,opts);

mu = qhat(1);
sigma = exp(qhat(2));
alpha = qhat(3);
nu = exp(qhat(4));

% Common histogram bins and observed PDF
edges = linspace(-10,10,501);
[counts,edges] = histcounts(data,edges);
bw = edges(2)-edges(1);
centers = edges(1:end-1)+bw/2;
obs = counts/(numel(data)*bw);

% Evaluate candidate PDFs at bin centers
% Approximation is good for sufficiently narrow bins
x = centers;
names = {'Gaussian','Student t','Skew-normal','Skew-t'};
pdfs = [pdf(pdNormal,x(:)), ...
        pdf(pdT,x(:)), ...
        skewpdf(x(:)',p)', ...
        skewt_pdf(x(:)',mu,sigma,alpha,nu)'];

if all(data > 0)
    names{end+1} = 'Lognormal';
    pdfs(:,end+1) = pdf(pdLognormal,x(:));
end

% Histogram and fits
figure;
histogram(data,edges,'Normalization','pdf');
hold on;
cols = lines(size(pdfs,2));
for j = 1:size(pdfs,2)
    plot(x,pdfs(:,j),'LineWidth',1.6,'Color',cols(j,:));
end
xlabel('Pixel intensity');
ylabel('Probability density');
title('Distribution Fits');
legend(['Histogram',names],'Location','best');
grid on;
xlim([-10 10]);

% Residuals: observed PDF minus fitted PDF
res = obs(:)-pdfs;

figure;
tiledlayout(size(pdfs,2),1,'TileSpacing','compact');
for j = 1:size(pdfs,2)
    nexttile;
    plot(x,res(:,j),'-','Color',cols(j,:),'LineWidth',1);
    hold on;
    yline(0,'k--');
    xlim([-10 10]);
    ylabel('Residual');
    title(names{j});
    grid on;
end
xlabel('Pixel intensity');

% -------- Local functions --------
function y = skewt_pdf(x,mu,sigma,alpha,nu)
    z = (x-mu)./sigma;
    arg = alpha.*z.*sqrt((nu+1)./(nu+z.^2));
    y = (2./sigma).*tpdf(z,nu).*tcdf(arg,nu+1);
end

function nll = skewt_nll(data,q)
    mu = q(1);
    sigma = exp(q(2));
    alpha = q(3);
    nu = exp(q(4));

    if ~isfinite(sigma) || ~isfinite(nu) || ...
       sigma < 1e-12 || nu < 0.05 || nu > 1e5 || ...
       abs(alpha) > 1e3
        nll = 1e100;
        return;
    end

    z = (data-mu)./sigma;
    arg = alpha.*z.*sqrt((nu+1)./(nu+z.^2));
    logp = log(2)-log(sigma) ...
        + log(max(tpdf(z,nu),realmin)) ...
        + log(max(tcdf(arg,nu+1),realmin));

    nll = -sum(logp);
    if ~isfinite(nll), nll = 1e100; end
end
