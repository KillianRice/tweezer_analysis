function funcOut = exponentialfit(analyVar, indivDataset, avgDataset)
    
    %% Fit_Template - Joe Whalen 2019.10.02
    
    % This function calls the base_fit script that does all of the fitting
    % and displays all of the plots for a given fitting routine. To use
    % this template, first save it as a new script with the name you want
    % to give to your fit function. Next, fill in your fit form with an
    % anonymous function handle with argumnents coeffs and x, coeffs is the
    % vector of free fit parameters and x is the vector of x coordinates
    % over which the function will be evaluated. Declare the independent
    % and dependent variables to be fitted with indVarField and
    % depVarField. These variables are strings that are the names of a
    % field in indivDataset. Typically the indVarField is imagevcoAtom, the
    % variable that was scanned during the experiment.

    % Gaussian with background and linear background
    
    form = @(coeffs, x) coeffs(1) .* exp(-2*(x-coeffs(2)).^2 ./ (coeffs(3)^2)) + coeffs(4) + coeffs(5).*x + coeffs(6) .* exp(-2*(x-coeffs(7)).^2 ./ (coeffs(8)^2)); % A * Exp[-(x-B)^2/(2 * C)^2] +  D + Ex
    
    indVarField = 'ImgPixelLengthList'; % independent variable
    depVarField1 = 'ImgXSum'; % dependent variable
    depVarField2 = 'ImgYSum'; % dependent variable
    coeffParam = 3;

    
    %% initial guess code
    % fill in this function to estimate the values of the fit parameters,
    % alternatively you can just have thist function return a constant
    % vector if you don't have a simple way of obtaining an initial guess
    function initialguess = x0(xdata, ydata)
        % code that guesses initial params
        % can also return a constant vector with length equal to the number
        % of parameters in the fit function
        
            
            %x0(1) = (min(y)-max(y)); %% For negative gaussians
            initialguess(1) = max(ydata);             %% For positive
            initialguess(2) =  sum(xdata.*ydata)/sum(ydata);
            %initialguess(2) = 81.96;
            initialguess(3) = sqrt(sum((xdata-initialguess(2)).^2.*ydata)/sum(ydata));
            initialguess(3) = 2;
            %x0(4) = max(y);            %% For negative gaussians
            initialguess(4) = min(ydata);             %% For positive
            initialguess(5) = 0;
            initialguess(6) = max(ydata)/2;
            initialguess(7) = sum(xdata.*ydata)/sum(ydata);
            initialguess(8) = initialguess(3)*10;
    end

    function initialguess = y0(xdata, ydata)
        % code that guesses initial params
        % can also return a constant vector with length equal to the number
        % of parameters in the fit function
        
            
            %x0(1) = (min(y)-max(y)); %% For negative gaussians
            initialguess(1) = max(ydata);             %% For positive
            initialguess(2) =  sum(xdata.*ydata)/sum(ydata);
            %initialguess(2) = 81.96;
            initialguess(3) = sqrt(sum((xdata-initialguess(2)).^2.*ydata)/sum(ydata));
            initialguess(3) = 2;
            %x0(4) = max(y);            %% For negative gaussians
            initialguess(4) = min(ydata);             %% For positive
            initialguess(5) = 0;
            initialguess(6) = max(ydata)/2;
            initialguess(7) = sum(xdata.*ydata)/sum(ydata);
            initialguess(8) = initialguess(3)*5;
    end


    %% options
    % the base_fit function is designed to be flexible, and can accept a
    % lot of different parameters to adjust 
    %%%%% Default Options %%%%%
        % options = struct(...
        % 'DataPlotFunction', @defaultDataPlot,...
        % 'AvgDataPlotFunction', @defaultAvgDataPlot,...
        % 'FitLinePlotFunction', @defaultFitLinePlot,...
        % 'AnnotateFunction', @defaultAnnotate,...
        % 'IndivFitPlotFunction', @defaultIndivFitPlot,...
        % 'PlotIndivFits' , false,...
        % 'PlotAvgFits', true,...
        % 'XAxisLabel', indVarField ,...
        % 'YAxisLabel', depVarField,...
        % 'FitLB', [],...
        % 'FitUB', [],...
        % 'FitOptions', struct('Display','off'),...
        % 'PlotInitialGuess', true, ...
        % 'InitialGuessPlotFunction', @defaultInitialGuessPlot,...
        % 'PlotAll', true,...
        % 'PlotAllAvgs', true, ...
        % 'CoeffNames', {{}},...
        % 'CoeffUnits', {{}},...
        % 'YAxisScale', 'linear',...
        % 'XAxisScale', 'linear'),...
        % %'FitTitle', (call from dbstack),...
        % 'Statistics', 'gaussian');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % THIS FUNCTION DOES NOT WORK USING INDIVIDUAL FITS
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAvgFits', true,...
        'PlotAll', false,...
        'PlotAllAvgs', true,...
        'PlotInitialGuess', true,...
        'Statistics', 'gaussian', ...
        'XAxisScale', 'linear');
    

    [xav1,yav1,yer1,coefflist1,coefflist_err1,avg_fit_coeffs_twzr1] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField1, @x0, options);
    [xav2,yav2,yer2,coefflist2,coefflist_err2,avg_fit_coeffs_twzr2] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField2, @y0, options);

    %% ============================================================
    %  EXTRACT FITTED PARAMETERS
    % ============================================================
    scanIDs = analyVar.uniqScanList;
    coeff2plotX = zeros(length(scanIDs), 1);
    coeff2plotY = zeros(length(scanIDs), 1);
    coeff2plotXerr = zeros(length(scanIDs), 1);
    coeff2plotYerr = zeros(length(scanIDs), 1);
    coeff2plotMeanAmp = zeros(length(scanIDs), 1);
    ratioSmallToLargeX = zeros(length(scanIDs), 1);
    ratioSmallToLargeY = zeros(length(scanIDs), 1);
    for id = 1:length(analyVar.uniqScanList)

        %check which gaussian is the smaller one
        if coefflist1{id}(3) < coefflist1{id}(8)
            smallerWidthSigIdx = 3;
            smallerWidthAmpIdx = 1;
            largerWidthSigIdx = 8;
            largerWidthAmpIdx = 6;
        else
            smallerWidthSigIdx = 8;
            smallerWidthAmpIdx = 6;
            largerWidthSigIdx = 3;
            largerWidthAmpIdx = 1;
        end

        if coefflist2{id}(3) < coefflist2{id}(8)
            smallerWidthSigIdy = 3;
            smallerWidthAmpIdy = 1;
            largerWidthSigIdy = 8;
            largerWidthAmpIdy = 6;
        else
            smallerWidthSigIdy = 8;
            smallerWidthAmpIdy = 6;
            largerWidthSigIdy = 3;
            largerWidthAmpIdy = 1;
        end

        fittedSigX = coefflist1{id}(smallerWidthSigIdx);
        fittedSigY = coefflist2{id}(smallerWidthSigIdy);
        fittedSigXLarge = coefflist1{id}(largerWidthSigIdx);
        fittedSigYLarge = coefflist2{id}(largerWidthSigIdy);

        AmpX = coefflist1{id}(smallerWidthAmpIdx) ./ ...
                (sqrt(2*pi) .* fittedSigY);
        AmpY = coefflist2{id}(smallerWidthAmpIdy) ./ ...
                (sqrt(2*pi) .* fittedSigX);
        AmpX_err = AmpX .* sqrt((coefflist_err1{id}(smallerWidthSigIdx) ./ fittedSigX).^2 + (coefflist_err1{id}(smallerWidthAmpIdx)./ AmpX).^2);
        AmpY_err = AmpY .* sqrt((coefflist_err2{id}(smallerWidthSigIdy) ./ fittedSigY).^2 + (coefflist_err2{id}(smallerWidthAmpIdy)./ AmpY).^2);

        AmpXLarge = coefflist1{id}(largerWidthAmpIdx);
        AmpYLarge = coefflist2{id}(largerWidthAmpIdy);       
        
        % Amplitude conversion else just store normal coefficient
        if coeffParam == 1
            coeff2plotX(id) = AmpX;
            coeff2plotY(id) = AmpY;
            coeff2plotXerr(id) = AmpX_err;
            coeff2plotYerr(id) = AmpY_err;
        else
            coeff2plotX(id) = analyVar.sizefactor * coefflist1{id}(smallerWidthSigIdx) * 1000000;
            coeff2plotY(id) = analyVar.sizefactor * coefflist2{id}(smallerWidthSigIdy) * 1000000;
            coeff2plotXerr(id) = analyVar.sizefactor * coefflist_err1{id}(smallerWidthSigIdx) * 1000000;
            coeff2plotYerr(id) = analyVar.sizefactor * coefflist_err2{id}(smallerWidthSigIdy) * 1000000;
        end

        coeff2plotMeanAmp(id) = mean([AmpX AmpY]);
        coeff2plotMeanAmpErr(id) = sqrt(AmpX_err.^2 + AmpY_err.^2)/2;
        
        
        ratioSmallToLargeX(id) = (coefflist1{id}(smallerWidthAmpIdx)*fittedSigX)/(AmpXLarge*fittedSigXLarge);
        ratioXerr(id) = abs(ratioSmallToLargeX(id)) * sqrt( ...
                        (coefflist_err1{id}(smallerWidthAmpIdx) / coefflist1{id}(smallerWidthAmpIdx))^2 + ...
                        (coefflist_err1{id}(smallerWidthSigIdx) / fittedSigX)^2 + ...
                        (coefflist_err1{id}(largerWidthAmpIdx) / AmpXLarge)^2 + ...
                        (coefflist_err1{id}(largerWidthSigIdx) / fittedSigXLarge)^2 );
        ratioSmallToLargeY(id) = (coefflist2{id}(smallerWidthAmpIdy)*fittedSigY)/(AmpYLarge*fittedSigYLarge);
        ratioYerr(id) = abs(ratioSmallToLargeY(id)) * sqrt( ...
                        (coefflist_err2{id}(smallerWidthAmpIdy)/ coefflist2{id}(smallerWidthAmpIdy))^2 + ...
                        (coefflist_err2{id}(smallerWidthSigIdy) / fittedSigY)^2 + ...
                        (coefflist_err2{id}(largerWidthAmpIdy) / AmpYLarge)^2 + ...
                        (coefflist_err2{id}(largerWidthSigIdy) / fittedSigYLarge)^2 );
    end

    %% ==============================================
    % Plotting
    figure;
    hold on;
    errorbar(scanIDs, coeff2plotX, coeff2plotXerr, 'o', 'LineWidth', 1.5,'MarkerFaceColor','b','MarkerEdgeColor','b');
    errorbar(scanIDs, coeff2plotY, coeff2plotYerr, 's', 'LineWidth', 1.5,'MarkerFaceColor','r','MarkerEdgeColor','r');
    ylabel(sprintf('Radius of Gaussian (um)'));
    hold off;
    
    % Labels and Legend
    xlabel('Scan ID');
    legend('Xsum', 'Ysum', 'Location', 'NorthWest');
    title(sprintf('Effective Pixelsize %d um/px',analyVar.sizefactor * 1000000))
    grid on
    %% ==============================================    
     % Plotting
    figure;
    hold on;
    errorbar(scanIDs, coeff2plotMeanAmp, coeff2plotMeanAmpErr, 'v', 'LineWidth', 1.5,'MarkerFaceColor','g','MarkerEdgeColor','g');
    ylabel(sprintf('1D Image Reduction Amplitude fit'));
    hold off;
    
    % Labels and Legend
    xlabel('Scan ID');
    legend('Fitted Amplitdue', 'Location', 'NorthWest');
    grid on

    %% ==============================================
    figure;
    hold on;

    errorbar(scanIDs, ratioSmallToLargeX, ratioXerr, 'o', 'LineWidth', 1.5,'MarkerFaceColor','b');
    errorbar(scanIDs, ratioSmallToLargeY, ratioYerr, 's', 'LineWidth', 1.5,'MarkerFaceColor','r');
    ylabel(sprintf('1D Image Reduction Ratio Small:Large Gaussian Area'));
    
    % Labels and Legend
    xlabel('Scan ID');
    legend('Xsum', 'Ysum', 'Location', 'NorthWest');
    grid on
    

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end

