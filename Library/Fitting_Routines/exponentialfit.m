function funcOut = exponentialfit(analyVar, indivDataset, avgDataset)
    
    %% Fit_Template - Joe Whalen 2019.10.02
    % Updated on 2025.02.04 by SKK to modify the form to A * Exp[-t/tau] + C 
    % Updated on 2026.05.21 by BLT to include normalization
    
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
    
    form = @(coeffs, x) coeffs(1)*exp(-x/coeffs(2)) + coeffs(3); % A * Exp[-t/tau] + C 
    
    indVarField = 'imagevcoAtom'; % independent variable
    %depVarField = 'sfiIntegral'; % dependent variable
    %depVarField = 'numberAtom'; % dependent variable
    depVarField = 'OD_TotalCountsImg1Raw';

    %normalize option (also adjust within initial guess)
    normalize = 0;
    
    %% initial guess code
    % fill in this function to estimate the values of the fit parameters,
    % alternatively you can just have thist function return a constant
    % vector if you don't have a simple way of obtaining an initial guess
    function initialguess = x0(xdata, ydata)
        % code that guesses initial params
        % can also return a constant vector with length equal to the number
        % of parameters in the fit function
        
        %normalize option
        normal = 0;
        if normal
            initialguess(1) = max(ydata)/max(ydata);
            initialguess(2) = (max(xdata)-min(xdata));
            initialguess(3) = min(ydata)/max(ydata);
        else
            initialguess(1) = 5000;
            initialguess(2) = (max(xdata)-min(xdata));
            initialguess(3) = min(ydata);
        end
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
    
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', true,...
        'PlotInitialGuess', true,...
        'Statistics', 'gaussian', ...
        'XAxisScale', 'linear');
    
    if normalize
        normalized_base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @x0, options)
    else
        base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @x0, options)
    end

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end

