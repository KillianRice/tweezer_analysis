function [xav,yav,yer,coefflist,coefflist_err] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, x0, useroptions)
    
    st = dbstack;
    call = st(2).name;
    %%%%% set up options and plotting functions %%%%%
    options = struct(...
        'DataPlotFunction', @defaultDataPlot,...
        'AvgDataPlotFunction', @defaultAvgDataPlot,...
        'FitLinePlotFunction', @defaultFitLinePlot,...
        'AnnotateFunction', @defaultAnnotate,...
        'IndivFitPlotFunction', @defaultIndivFitPlot,...
        'PlotIndivFits' , true,...
        'PlotAvgFits', true,...
        'XAxisLabel', indVarField ,...
        'YAxisLabel', depVarField,...
        'FitLB', [],...
        'FitUB', [],...
        'FitOptions', struct('Display','off'),...
        'PlotInitialGuess', true, ...
        'InitialGuessPlotFunction', @defaultInitialGuessPlot,...
        'PlotAll', true,...
        'PlotAllAvgs', true, ...
        'CoeffNames', {{}},...
        'CoeffUnits', {{}},...
        'YAxisScale', 'linear',...
        'XAxisScale', 'linear',...
        'FitTitle', call,...
        'Statistics', 'gaussian');
    
    if nargin > 7
        opts = fieldnames(useroptions);
        for i = 1:numel(opts)
            options.(opts{i}) = useroptions.(opts{i});
        end
    end
    
    myDataPlot = options.('DataPlotFunction');
    myAvgDataPlot = options.('AvgDataPlotFunction');
    myFitLinePlot = options.('FitLinePlotFunction');
    myAnnotate = options.('AnnotateFunction');
    myInitialPlot = options.('InitialGuessPlotFunction');
    
    plotIndivFits = options.PlotIndivFits;
    plotAvgFits = options.PlotAvgFits;
    plotInitialGuess = options.PlotInitialGuess;
    plotAll = options.PlotAll;
    plotAllAvgs = options.PlotAllAvgs;

    fitLB = options.FitLB;
    fitUB = options.FitUB;
    fitOptions = options.FitOptions;
    
    xlabeltext = options.XAxisLabel;
    ylabeltext = options.YAxisLabel;
    
    coeffNames = options.CoeffNames;
    coeffUnits = options.CoeffUnits;
    
    xAxisScale = options.XAxisScale;
    yAxisScale = options.YAxisScale;
    
    fitTitle = options.FitTitle;
    
    weighting = options.Statistics;
    
    
    if plotIndivFits
    
        [xdata, ydata] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
        coeffs = cell(analyVar.numBasenamesAtom,1);
        uncs = cell(analyVar.numBasenamesAtom,1);

        for i = 1:analyVar.numBasenamesAtom

            % try to correct for data that are not the same size
            if size(xdata{i}) ~= size(ydata{i})
                warning(['Dimensions of xdata, ydata not the same. ' ...
                    'Trying to fix, but may lead to unpredictable results.'])
                ydata{i} = ydata{i}';
            end

            % fit the data
            initialguess = x0(xdata{i},ydata{i});
            [coeffs{i},~,~,CovB,rchisq,~] = nlinfit(xdata{i},ydata{i},form,initialguess);
            uncs{i} = sqrt(diag(CovB)); % 1 sigma uncertainty from covariance matrix
            
            % plot the data
            fitx = linspace(min(xdata{i}),max(xdata{i}),1000);
            
            figure
            hold on
                if plotInitialGuess
                    defaultInitialGuessPlot(fitx, form(initialguess, fitx), i, analyVar);
                end
                myDataPlot(xdata{i},ydata{i},i,analyVar);
                myFitLinePlot(fitx, form(coeffs{i},fitx),i,analyVar);
                myAnnotate(coeffs{i},uncs{i}, coeffNames, coeffUnits);
                disp(strcat(['Fit data for ' num2str(analyVar.timevectorAtom(i))]))
                disp(coeffs{i})
                xlabel(xlabeltext,'Interpreter','none');
                ylabel(ylabeltext,'Interpreter','none');
                legend(num2str(analyVar.timevectorAtom(i)));
                set(gca, 'YScale', yAxisScale);
                set(gca, 'XScale', xAxisScale);
                title(strcat([fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' \nu = ',...
                    num2str(length(ydata{i})-length(coeffs{i}))]));
            hold off
        end
        
        if plotAll
            figure
            hold on
            for i = 1:analyVar.numBasenamesAtom
                myDataPlot(xdata{i},ydata{i},i,analyVar);
                %myFitLinePlot(fitx, form(coeffs{i},fitx),i,analyVar);
                xlabel(xlabeltext);
                ylabel(ylabeltext);
            end
            legend(num2str(analyVar.timevectorAtom));
            set(gca, 'YScale', yAxisScale);
            set(gca, 'XScale', xAxisScale);
            hold off
        end
        
    end
    
    if length(analyVar.timevectorAtom) > 1 && plotAvgFits
        
        [xavg, yavg, yerr] = get_averages(analyVar, indivDataset, avgDataset,...
            indVarField, depVarField, weighting);
        scanIDs = analyVar.uniqScanList;
        avg_coeffs = cell(length(scanIDs),1);
        avg_unc = cell(length(scanIDs),1);
                
        for i = 1:length(scanIDs)
            
            if size(xavg{i}) ~= size(yavg{i})
                warning(['Dimensions of xdata, ydata not the same. ' ...
                    'Trying to fix, but may lead to unpredictable results.'])
                yavg{i} = yavg{i}';
            end
            
            % fit the data
            initialguess = x0(xavg{i}, yavg{i});
            
            weights = 1./(yerr{i} + 1).^2;
            [avg_coeffs{i},~,~,CovB,rchisq,~] = nlinfit(xavg{i},yavg{i},form,initialguess,...
                'Weights',weights);
            avg_unc{i} = sqrt(diag(CovB));
            % plot the data
            fitx = linspace(min(xavg{i}),max(xavg{i}),1000);
            
            figure
            hold on
                if plotInitialGuess
                    defaultInitialGuessPlot(fitx, form(initialguess, fitx), i, analyVar);
                end
                myAvgDataPlot(xavg{i},yavg{i},yerr{i},i,analyVar);
                myFitLinePlot(fitx, form(avg_coeffs{i},fitx),i,analyVar);
                myAnnotate(avg_coeffs{i}, avg_unc{i}, coeffNames, coeffUnits);
                xlabel(xlabeltext,'Interpreter','none');
                ylabel(ylabeltext,'Interpreter','none');
                legend(num2str(scanIDs(i)));
                set(gca, 'YScale', yAxisScale);
                set(gca, 'XScale', xAxisScale);
                title(strcat([fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' \nu = ',...
                    num2str(length(xavg{i})-length(avg_coeffs{i}))]));
            hold off
            
            
        end
        
        if plotAllAvgs
            figure
            hold on
            for i = 1:length(scanIDs)
                myAvgDataPlot(xavg{i},yavg{i},yerr{i},i,analyVar);
                myFitLinePlot(fitx, form(avg_coeffs{i},fitx),i,analyVar);
                xlabel(xlabeltext,'Interpreter','none');
                ylabel(ylabeltext,'Interpreter','none');
                legend(num2str(scanIDs(i)));
                set(gca, 'YScale', yAxisScale);
                set(gca, 'XScale', xAxisScale);
            end
            legend(num2str(scanIDs));
            hold off
        end
    xav = xavg;
    yav = yavg;
    yer = yerr;
    coefflist = avg_coeffs;    
    coefflist_err = avg_unc;
    end
  
end
    
function h = defaultDataPlot(x,y,i,analyVar)
    h = plot(x,y,...
    'LineStyle','none',...
    'Marker', 'o',...
    'MarkerSize', analyVar.markerSize,...
    'MarkerFaceColor', analyVar.COLORS(i,:),...
    'MarkerEdgeColor', 'k',...
    'Color', analyVar.COLORS(i,:));
end

function h = defaultAvgDataPlot(x,y,yerr,i,analyVar)
    h = errorbar(x,y,yerr,...
        'LineStyle','none',...
        'Marker', 'o',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(i,:),...
        'MarkerEdgeColor', 'k',...
        'Color', analyVar.COLORS(i,:));
end

function h = defaultFitLinePlot(x,y,i,analyVar)
    h = plot(x,y,...
        'LineStyle', '-',...
        'Marker', 'none',...
        'LineWidth', 1,...
        'Color', analyVar.COLORS(i,:));
end

function an = defaultAnnotate(coeffs, err, coeffNames, coeffUnits)
    
    dim = [.7 .5 .3 .3];
    
    if isempty(coeffNames)
        for i = 1:numel(coeffs)
            coeffNames{i} = ['Coeff ', num2str(i)];
        end
    end
    
    if isempty(coeffUnits)
        for i = 1:numel(coeffs)
            coeffUnits{i} = '';
        end
    end
    
    strs = cell(numel(coeffs),1);
    for i = 1:numel(coeffs)
        if i < numel(coeffs)
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}, newline];
        else
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}];
        end
    end
    
    an = annotation('textbox', dim, 'String', strjoin(strs),...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
end

function h = defaultInitialGuessPlot(x,y,i,analyVar)
    h = plot(x,y,...
    'LineStyle','--',...
    'Color', analyVar.COLORS(i+2,:));
end

