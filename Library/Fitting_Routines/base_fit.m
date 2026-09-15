function [xav,yav,yer,coefflist,coefflist_err,avg_fit_coeffs_twzr] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, x0, useroptions)
    
    %%default outputs
    xav = 0;
    yav = 0;
    yer = 0;
    coefflist = 0;    
    coefflist_err = 0;

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
        'Statistics', 'gaussian', ...
        'PlotNormalized', false, ...
        'PlotNormToMeasuredParams', false);
    
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
    plotNormalized = options.PlotNormalized;
    plotNormToMeasuredParams = options.PlotNormToMeasuredParams; %% Include For normalizing fits and plots to params: Atom Num, Power, etc...

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
     
    %% Cycle through all tweezer spots, if not using tweezers then just do once
    if analyVar.UseTweezer
        load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
        numTweezers = size(tweezerROI.centersXY,1);
    else
        numTweezers = 1;
    end

     %% Creating stuct to save COEFFICIENTS of fits for later plotting an analysis
     avg_fit_coeffs_twzr = struct;
     avg_fit_coeffs_twzr.scanIDs = cell(numTweezers,1);
     avg_fit_coeffs_twzr.coeffs = cell(numTweezers,1);
     avg_fit_coeffs_twzr.uncs = cell(numTweezers,1);
     avg_fit_coeffs_twzr.numTweezers = numTweezers;


     %% ============================================================================
     %              BEGIN FITTING
     %  ----------------------------------------------------------------------------

    for tweezerNum = 1:numTweezers
        if analyVar.UseTweezer
            fprintf('plotting tweezers %d', tweezerNum);
        end

        if analyVar.SkipIndivTwzrFitting && analyVar.UseTweezer
            fprintf('Skipped tweezers %d', tweezerNum);
            continue;
        end

        %% ============================================================================
        %               FITS FOR EACH INDIVDATASET
        %  ----------------------------------------------------------------------------
        if plotIndivFits
        
            [xdata, ydata] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset, tweezerNum);   %%% Added Tweezer capability
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
                    if plotNormalized                   %% Edit for plotting normalized data and fits
                        if plotInitialGuess
                            defaultInitialGuessPlot(fitx, form(initialguess, fitx)/initialguess(end), i, analyVar);
                        end
                        myDataPlot(xdata{i},ydata{i}/coeffs{i}(end),i,analyVar);
                        myFitLinePlot(fitx, form(coeffs{i},fitx)/coeffs{i}(end),i,analyVar);
                    else
                        if plotInitialGuess
                            defaultInitialGuessPlot(fitx, form(initialguess, fitx), i, analyVar);
                        end
                        myDataPlot(xdata{i},ydata{i},i,analyVar);
                        myFitLinePlot(fitx, form(coeffs{i},fitx),i,analyVar);
                    end
                    myAnnotate(coeffs{i},uncs{i}, coeffNames, coeffUnits);
                    % disp(strcat(['Fit data for ' num2str(analyVar.timevectorAtom(i))]))
                    % disp(coeffs{i})
                    xlabel(xlabeltext,'Interpreter','none');
                    ylabel(ylabeltext,'Interpreter','none');
                    legend(num2str(analyVar.timevectorAtom(i)));
                    set(gca, 'YScale', yAxisScale);
                    set(gca, 'XScale', xAxisScale);
                    if analyVar.UseTweezer
                        title(strcat([fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' A\nu = ',...
                            num2str(length(ydata{i})-length(coeffs{i})), ' Tweezer ROI: ', tweezerNum]));
                    else
                        title(strcat([fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' \nu = ',...
                            num2str(length(ydata{i})-length(coeffs{i}))]));
                    end
                hold off
            end
        
            %% ============================================================
            %          PLOT ALL INDIVDATASET FITS ONTO SINGLE PLOT
            if plotAll
                figure
                hold on
                for i = 1:analyVar.numBasenamesAtom
                    if plotNormalized
                        myDataPlot(xdata{i},ydata{i}/coeffs{i}(end),i,analyVar);
                    else
                        myDataPlot(xdata{i},ydata{i},i,analyVar);
                    end
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
        
        %% ============================================================================
        %         FITS FOR DATA AVERAGED BASED ON SAME IMAGEVCOATOM WITHIN
        %                   INDIVDATASETS WITH SAME SCAN IDS
        %  ----------------------------------------------------------------------------
        
        if length(analyVar.timevectorAtom) > 1 && plotAvgFits
        
            [xavg, yavg, yerr] = get_averages(analyVar, indivDataset, avgDataset,...
                indVarField, depVarField, weighting, tweezerNum);                   %%%% Change for including separate tweezers
            scanIDs = analyVar.uniqScanList;
            avg_coeffs = cell(length(scanIDs),1);
            avg_unc = cell(length(scanIDs),1);
            


            %% Normalize data based on Params (E.g. AtomNum, Power, etc..)
            if plotNormToMeasuredParams
                [num, ~, tx, ~, ty, ~] = get_num_temp_averages(analyVar, indivDataset);
                [spec413,~,spec461,~] = get_daq_averages(analyVar, indivDataset);
                disp(['spec 413: ', num2str(spec413)]);
                disp(['spec 461: ', num2str(spec461)]);
                disp(['num: ', num2str(num)]);
            end
        
            for i = 1:length(scanIDs)
        
                if size(xavg{i}) ~= size(yavg{i})
                    warning(['Dimensions of xdata, ydata not the same. ' ...
                        'Trying to fix, but may lead to unpredictable results.'])
                    yavg{i} = yavg{i}';
                end

                if plotNormToMeasuredParams
                    yavg{i} = yavg{i}/num(i)/spec413(i)/spec461(i);
                    yerr{i} = yerr{i}/num(i)/spec413(i)/spec461(i);
                end
        
                % fit the data
                initialguess = x0(xavg{i}, yavg{i});
        
                weights = 1./(yerr{i} + 1).^2;
                %disp(yerr{i})
                [avg_coeffs{i},~,~,CovB,rchisq,~] = nlinfit(xavg{i},yavg{i},form,initialguess,...
                    'Weights',weights);
                avg_unc{i} = sqrt(diag(CovB));
                % plot the data
                fitx = linspace(min(xavg{i}),max(xavg{i}),1000);
        
                figure
                hold on
                    if plotNormalized                   %% Edit for plotting normalized data and fits
                        if plotInitialGuess
                            defaultInitialGuessPlot(fitx, form(initialguess, fitx)/initialguess(end), i, analyVar);
                        end
                        % disp(avg_coeffs{i}(end))
                        myAvgDataPlot(xavg{i},yavg{i}/avg_coeffs{i}(end),yerr{i}/avg_coeffs{i}(end),i,analyVar);
                        myFitLinePlot(fitx, form(avg_coeffs{i},fitx)/avg_coeffs{i}(end),i,analyVar);
                    else
                        if plotInitialGuess
                            defaultInitialGuessPlot(fitx, form(initialguess, fitx), i, analyVar);
                        end
                        myAvgDataPlot(xavg{i},yavg{i},yerr{i},i,analyVar);
                        myFitLinePlot(fitx, form(avg_coeffs{i},fitx),i,analyVar);
                    end
                    myAnnotate(avg_coeffs{i}, avg_unc{i}, coeffNames, coeffUnits);
                    xlabel(xlabeltext,'Interpreter','none');
                    ylabel(ylabeltext,'Interpreter','none');
                    legend(num2str(scanIDs(i)),'Data','Fit');
                    set(gca, 'YScale', yAxisScale);
                    set(gca, 'XScale', xAxisScale);
                    if analyVar.UseTweezer
                        title(strcat([' Tweezer: ', num2str(tweezerNum), ' ', fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' A\nu = ',...
                            num2str(length(yavg{i})-length(avg_coeffs{i}))]));
                    else
                        title(strcat([fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' \nu = ',...
                            num2str(length(yavg{i})-length(avg_coeffs{i}))]));
                    end
                hold off
        
        
            end

            %Creating stuct to save fitted coefficents for later plotting
            avg_fit_coeffs_twzr.scanIDs{tweezerNum} = scanIDs;
            avg_fit_coeffs_twzr.coeffs{tweezerNum} = avg_coeffs;
            avg_fit_coeffs_twzr.uncs{tweezerNum} = avg_unc;
        
            %% ============================================================
            %          PLOT ALL AVERAGED FITS ONTO SINGLE PLOT
            if plotAllAvgs
                % figure
                % hold on
                % for i = 1:length(scanIDs)
                %     myAvgDataPlot(xavg{i},yavg{i},yerr{i},i,analyVar);
                %     myFitLinePlot(fitx, form(avg_coeffs{i},fitx),i,analyVar);
                %     xlabel(xlabeltext,'Interpreter','none');
                %     ylabel(ylabeltext,'Interpreter','none');
                %     legend(num2str(scanIDs(i)));
                %     set(gca, 'YScale', yAxisScale);
                %     set(gca, 'XScale', xAxisScale);
                % end
                % legend(num2str(scanIDs));
                % hold off
                figure;
                hold on;
                
                legendHandles = gobjects(length(scanIDs),1);
                legendLabels  = cell(length(scanIDs),1);
                
                for i = 1:length(scanIDs)
                
                    if plotNormalized
                        % Plot data and keep its graphics handle for the legend
                        legendHandles(i) = myAvgDataPlot( ...
                            xavg{i}, ...
                            yavg{i}/avg_coeffs{i}(end), ...
                            yerr{i}/avg_coeffs{i}(end), ...
                            i, ...
                            analyVar);
                    
                        % Plot the corresponding fit, but do not add it separately to the legend
                        fitHandle = myFitLinePlot( ...
                            fitx, ...
                            form(avg_coeffs{i},fitx)/avg_coeffs{i}(end), ...
                            i, ...
                            analyVar);
                    else
                        % Plot data and keep its graphics handle for the legend
                        legendHandles(i) = myAvgDataPlot( ...
                            xavg{i}, ...
                            yavg{i}, ...
                            yerr{i}, ...
                            i, ...
                            analyVar);
                    
                        % Plot the corresponding fit, but do not add it separately to the legend
                        fitHandle = myFitLinePlot( ...
                            fitx, ...
                            form(avg_coeffs{i},fitx), ...
                            i, ...
                            analyVar);
                    end
                
                    fitHandle.HandleVisibility = 'off';
                
                    legendLabels{i} = num2str(scanIDs(i));
                end
                
                xlabel(xlabeltext,'Interpreter','none');
                ylabel(ylabeltext,'Interpreter','none');
                
                set(gca,'YScale',yAxisScale);
                set(gca,'XScale',xAxisScale);
                
                legend( ...
                    legendHandles, ...
                    legendLabels, ...
                    'Location','best');
                
                hold off;
            end
        xav = xavg;
        yav = yavg;
        yer = yerr;
        coefflist = avg_coeffs;    
        coefflist_err = avg_unc;
        end
    end

    %% =====================================================================
    %
    %     DO A BASEFIT FITTING FOR THE AVERAGE OF ALL TWEEZER SPOTS
    %         (ABOVE WAS DOWN FOR EACH INDIVIDUAL TWEEZER SPOT)
    if analyVar.UseTweezer && analyVar.avgTweezerBasefit
            scanIDs = analyVar.uniqScanList;
            coeffs = cell(length(scanIDs),1);
            uncs = cell(length(scanIDs),1);
            xdata = avgDataset.ODCountsX;
            ydata = avgDataset.ODCountsY;
            yerr = avgDataset.ODCountsErr;
    
            for i = 1:length(scanIDs)

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
                    errorbar( ...
                            xdata{i}, ...
                            ydata{i}, ...
                            yerr{i}, ...
                            'LineStyle','none',...
                            'Marker', 'o',...
                            'MarkerSize', analyVar.markerSize,...
                            'MarkerFaceColor', analyVar.COLORS(i,:),...
                            'MarkerEdgeColor', 'k',...
                            'Color', analyVar.COLORS(i,:));
                    myFitLinePlot(fitx, form(coeffs{i},fitx),i,analyVar);
                    myAnnotate(coeffs{i},uncs{i}, coeffNames, coeffUnits);
                    % disp(strcat(['Fit data for ' num2str(analyVar.timevectorAtom(i))]))
                    % disp(coeffs{i})
                    xlabel(xlabeltext,'Interpreter','none');
                    ylabel(ylabeltext,'Interpreter','none');
                    legend(num2str(analyVar.timevectorAtom(i)));
                    set(gca, 'YScale', yAxisScale);
                    set(gca, 'XScale', xAxisScale);
                    title(strcat(['Avg All Twzr', fitTitle, ' \chi^2_{\nu} = ',num2str(rchisq),' A\nu = ',...
                        num2str(length(ydata{i})-length(coeffs{i}))]));
                hold off
            end

                figure
                hold on
                        for i = 1:length(scanIDs)

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
 
                                                % Plot data and keep its graphics handle for the legend
                        legendHandles(i) = myAvgDataPlot( ...
                            xdata{i}, ...
                            ydata{i}, ...
                            yerr{i}, ...
                            i, ...
                            analyVar);
                    myFitLinePlot(fitx, form(coeffs{i},fitx),i,analyVar);
                    fitHandle.HandleVisibility = 'off';
                
                    legendLabels{i} = num2str(scanIDs(i));
                    %myAnnotate(coeffs{i},uncs{i}, coeffNames, coeffUnits);
             
                        end

                xlabel(xlabeltext,'Interpreter','none');
                ylabel(ylabeltext,'Interpreter','none');
                
                set(gca,'YScale',yAxisScale);
                set(gca,'XScale',xAxisScale);
                title(strcat(['Avg All Twzr', fitTitle]));
                
                legend( ...
                    legendHandles, ...
                    legendLabels, ...
                    'Location','best');
                        hold off
    end
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ==========================================
%       HELPER FUNCTIONS
%  ==========================================
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
    
    dim = [.32 .4 .3 .3];
    
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

