function funcOut = rabi_oscillation_fit(analyVar, indivDataset, avgDataset)
    
    %% lossspectrafit - Soumya K kanungo 2020.06.12
    % Fits MCS data to a typical loss profile.
    
    form = @(coeffs, x) coeffs(1) * ...
        (sin(coeffs(2).*x + coeffs(3))) .* exp(-(coeffs(4).*x))  + ...
        coeffs(5);% Decaying oscillation profile
    
    indVarField = 'imagevcoAtom'; % independent variable
    depVarField = 'sfiIntegral_roi1_ratio'; % dependent variable
    depVarField2 = 'sfiIntegral_roi2_ratio'; % dependent variable
    
    %% toggle on and off saving coefficients to an excel sheet (added by npi)
    saveVals = 0;

    %% toggle on and off contast calculation -- only for one Scan ID at a time for now (added by npi)
    calculate_contrast = 1;
    
    %% initial guess code
    function initialguess = x0(xdata, ydata)
        initialguess = [0.4,50*2*pi*1,0.5*pi,100,.4];
        %initialguess(1) = min(ydata)-max(ydata);
        initialguess(5) = mean(ydata);
    end
    

    %% compute the contrast based on the max and min of the averaged data
    if (calculate_contrast == 1)
        [xavg, yavg1, yerr1] = get_averages(analyVar, indivDataset, avgDataset,...
                indVarField, depVarField, 'gaussian');
        [xavg, yavg2, yerr2] = get_averages(analyVar, indivDataset, avgDataset,...
                indVarField, depVarField2, 'gaussian');
        %pad the arrays in case of array length mismatch
        % celldisp(yavg1)
        % if length(yavg1) < length(yavg2)
        %    yavg1(end+1 : length(yavg2)) = 0.5; 
        %     yerr1(end+1 : length(yerr2)) = 0.5;
        % elseif length(yavg2) < length(yavg1)
        %     yavg2(end+1 : length(yavg1)) = 0.5;
        %     yerr2(end+1 : length(yerr1)) = 0.5;
        % end
        % both_amps = {[yavg1{:}; yavg2{:}]};
        % both_errs = {[yerr1{:}; yerr2{:}]};

        %determine which oscillation roi is g.s. by taking the sum of the
        %points and seeing which one is greater
        if(sum(yavg1{:}) > sum(yavg2{:}))
            [max_amp, max_pos] = max(yavg1{:});
            [min_amp, min_pos] = min(yavg1{:});

            max_err = yerr1{1}(max_pos);
            min_err = yerr1{1}(min_pos);
            
        else
            [max_amp, max_pos] = max(yavg2{:});
            [min_amp, min_pos] = min(yavg2{:});

            max_err = yerr2{1}(max_pos);
            min_err = yerr2{1}(min_pos);

        end   
        contrast_val = (max_amp - min_amp)/(max_amp + min_amp);
        
        %compute error based on error propogation using the error of the two
        %points
        add_err = sqrt((max_err)^2 + (min_err)^2);
    
        contrast_err = sqrt((add_err/(max_amp + min_amp))^2 +...
            (add_err/(max_amp - min_amp))^2);

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
    
    %% If plotting Indiv Fits, contrast will no longer appear
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', false,...
        'XAxisLabel', 'time (ms)' ,...
        'YAxisLabel', 'MCS Signal Ratio',...
        'CoeffNames', {{'Ampl.', 'Freq.', 'phase', 'decay rate','Coh. time'}},...
        'CoeffUnits', {{'','MHz','rads','/micsec','micsec.'}},...
        'AnnotateFunction', @myAnnotate,...
        'PlotInitialGuess', true);
    [xav1,yav1,yer1,coefflist1,coefflist_err1]  = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @x0, options);
    if(calculate_contrast == 1)
        h = findall(gcf, 'Type', 'textboxshape');
        h.String{end} = ['Contrast: ', unc_string(contrast_val,contrast_err)];
    end 
    [xav2,yav2,yer2,coefflist2,coefflist_err2]  = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField2, @x0, options);
    if(calculate_contrast == 1)
        h = findall(gcf, 'Type', 'textboxshape');
        h.String{end} = ['Contrast: ', unc_string(contrast_val,contrast_err)];
    end 

    coefflist1new = coefflist2{1};
    target_y = 0.5;
    
    t_min = 0;
    t_max = 6e-3;   % choose your real physical max time here
    N = 10000;
    
    t_grid = linspace(t_min, t_max, N);
    
    root_fun = @(t) form(coefflist1new, t) - target_y;
    
    y_grid = arrayfun(root_fun, t_grid);
    
    % Find first crossing of target_y
    idx = find(y_grid(1:end-1).*y_grid(2:end) <= 0, 1, 'first');
    
    if isempty(idx)
        fprintf('No crossing found in the specified interval.\n');
    else
        bracket = [t_grid(idx), t_grid(idx+1)];
    
        earliest_time = fzero(root_fun, bracket);
    
        fprintf('Earliest solution: %.8g ms\n', earliest_time);
        fprintf('Earliest solution: %.4f ns\n', earliest_time*1e6);
    end
    
    if (saveVals== 1)
        filename = "MMWavePlotVals.xlsx";
        sheet    = "Sheet1";

        timeStr = strjoin(string(analyVar.timevectorAtom), ",");   % "1.2,3.4,5.6"
        plugInStr = strjoin(string(analyVar.plugInVec), ",");
        
        % Save several values to the array
        % time stamp | plug in vector | Amplitude | Frequency | Phase | Decay Rate | Coherence Time | Frequency/Decay | Contrast
        row = {timeStr, plugInStr, coefflist1{1}(1), coefflist1{1}(2)/10^3/2/pi, coefflist1{1}(3), coefflist1{1}(4)/10^3, 1/(coefflist1{1}(4)), coefflist1{1}(2)*coefflist1{1}(4), contrast_val, '58 s', '58 d', indivDataset{1,1}.fileAtom{1}};  % cell row
        
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

    dim = [.5 .7 .3 .3];
    coeffs(2) = coeffs(2)/10^3/2/pi  % converting the Omega to MHz units 
    err(2) = err(2)/10^3/2/pi      % converting the Omega to MHz units 
    strs = cell(6,1);
    coeffs(4) = coeffs(4)/10^3;  % converting rate to per microsecond units
    err(4) = err(4)/10^3;
    coeffs(5) = 1/(coeffs(4));
    for i = 1:5
        if i <= numel(coeffs)
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}, newline];
            strs{i}
        end
    end

    strs{6} = ['Freq/Decay: ', unc_string(coeffs(2)/coeffs(4), (coeffs(2)/coeffs(4))* (err(2)/coeffs(2) + err(4)/coeffs(4))), newline];
    
    
    an = annotation('textbox', dim, 'String', strjoin(strs),...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');

end

