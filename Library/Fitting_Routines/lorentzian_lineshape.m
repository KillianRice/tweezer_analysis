function funcOut = lorentzian_lineshape(analyVar, indivDataset, avgDataset)
    
    % Lorentzian fit - Soumya K kanungo 2020.07.31
    % Fits field to a typical lorentzian. 
    
    form = @(coeffs, x) coeffs(1)...
        .*(coeffs(3).^2)./((x-coeffs(2)).^2+coeffs(3).^2) + coeffs(4);

    %% toggle saveVals on and off to save values to an external excel sheet
    %% added by npi
    saveVals = 0 ;
    
    indVarField = 'imagevcoAtom'; % independent variable
    %depVarField = 'sfiIntegral'; % dependent variable
    depVarField = 'sfiIntegral_roi1_ratio';
    
    %% initial guess code
    function initialguess = x0(xdata, ydata)
        initialguess = [1,1048.9,0.1,.1];
        initialguess(1) = max(ydata)-min(ydata);
        initialguess(3) = (max(xdata)-min(xdata))/20;
        [maxvalue,position] = max(ydata);
        initialguess(2) = xdata(position);
        %initialguess(2) = 16663.5;
    end


    %% options
    % the base_fit function is designed to be flexible, and can accept a
    % lot of different parameters to adjust 
    %%%%% Default Options %%%%%
%         options = struct(...
%         'DataPlotFunction', @defaultDataPlot,...
%         'AvgDataPlotFunction', @defaultAvgDataPlot,...
%         'FitLinePlotFunction', @defaultFitLinePlot,...
%         'AnnotateFunction', @defaultAnnotate,...
%         'IndivFitPlotFunction', @defaultIndivFitPlot,...
%         'PlotIndivFits' , true,...
%         'PlotAvgFits', true,...
%         'XAxisLabel', indVarField ,...
%         'YAxisLabel', depVarField,...
%         'FitLB', [],...
%         'FitUB', [],...
%         'FitOptions', struct('Display','off'),...
%         'PlotInitialGuess', true, ...
%         'InitialGuessPlotFunction', @defaultInitialGuessPlot,...
%         'PlotAll', true,...
%         'PlotAllAvgs', true, ...
%         'CoeffNames', {{}},...
%         'CoeffUnits', {{}},...
%         'YAxisScale', 'linear',...
%         'XAxisScale', 'linear');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', false,...
        'XAxisLabel', 'WindFreak Synth (MHz)' ,...
        'CoeffNames', {{'Amp1.', '\mu_1.', '\sigma_1'}},...
        'CoeffUnits', {{'','MHz','MHz'}},...
        'AnnotateFunction', @myAnnotate,...
        'PlotInitialGuess', true);
    [xav,yav,yer,coefflist1,coefflist_err] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @x0, options)
    
    if (saveVals== 1)
        filename = "MMWavePlotVals.xlsx";
        sheet    = "Sheet1";

        timeStr = strjoin(string(analyVar.timevectorAtom), ",");   % "1.2,3.4,5.6"
        plugInStr = strjoin(string(analyVar.plugInVec), ",");
        
        % Save several values to the array
        % time stamp | plug in vector | Amplitude | Frequency | 
        row = {timeStr, plugInStr, coefflist1{1}(1), coefflist1{1}(2), coefflist1{1}(3), coefflist1{1}(4), '58 s', '58 d', indivDataset{1,1}.fileAtom{1}};  % cell row
        
        if isfile(filename)
            C = readcell(filename, "Sheet", sheet);
            nextRow = size(C,1) + 1;
        else
            nextRow = 1;
        end
   
        range = "A" + nextRow;
        writecell(row, filename, "Sheet", sheet, "Range", range);
    end 

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end
function an = myAnnotate(coeffs, err, coeffNames, coeffUnits)
    dim = [.7 .5 .3 .3];
    strs = cell(3,1);
    for i = 1:3
        if i <= numel(coeffs)
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}, newline];
            strs{i};
        end
    end
    
    an = annotation('textbox', dim, 'String', strjoin(strs),...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
end

