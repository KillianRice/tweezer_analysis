function funcOut = zeeman_sfi_triple_lorentzian(analyVar, indivDataset, avgDataset)
    
    %function that fits to a triple peaked Lorentzian
      %form = @(coeffs,x) coeffs(1) * exp(-(x-coeffs(6)-coeffs(2)).^2 ./ (2*coeffs(3)^2)) + coeffs(5) * exp(-(x-coeffs(6)).^2 ./ (2*coeffs(7)^2)) + coeffs(8) * exp(-(x-coeffs(6)+coeffs(2)).^2 ./ (2*coeffs(9)^2)) + coeffs(4);
    
    % coeffs = [A1, A0, A2, sigma1, sigma0, sigma2, Delta, x0, bg]
    %%NOTE the sqrt on the widths
    form = @(coeffs,x) ...
    coeffs(1)*(coeffs(4) ./ ((x-(coeffs(8)-coeffs(7))).^2 + coeffs(4)^2)) + ...
    coeffs(2)*(coeffs(5) ./ ((x-coeffs(8)).^2 + coeffs(5)^2)) + ...
    coeffs(3)*(coeffs(6) ./ ((x-(coeffs(8)+coeffs(7))).^2 + coeffs(6)^2)) + ...
    coeffs(9);

    indVarField = 'imagevcoAtom';
    depVarField = 'sfiIntegral_roi1_ratio';
    %depVarField = 'sfiIntegral';
    
    %% toggle saveVals on and off to save values to an external excel sheet
    %% added by npi
    saveVals = 0 ;
   
    function x0 = initial_guess(x, y)

        %Background
        bg = median(y(1:round(end*0.1))); 
        %y values minus background for other fitting
        y2 = y - bg;
        %Center Peak
        %x0 = sum(x.*y2)/sum(y2);
        %x0 = 21316.9;
        x0 = 20193.6;
        %other peak locations        %sort x values in ascending order for scans that are backwards
        [xs, idx] = sort(x);
        ys = y2(idx);
        [pks, locs] = findpeaks(ys, xs, ...
            'MinPeakProminence', 0.1*max(y2), ...
            'SortStr','descend');
        locs = sort(locs(1:3));
        %Delta = mean([x0 - locs(1), locs(3) - x0]);
        Delta = 3.6;
        disp(Delta)
        %Amplitudes
        A1 = interp1(x,y2,x0,'linear','extrap');
     
        A1 = 0.06;
        A0 = interp1(x,y2,x0-Delta,'linear','extrap');
        A0 = 0.035;
        A2 = interp1(x,y2,x0+Delta,'linear','extrap');
        A2 = 0.025;
        %Widths
        sigma_est = sqrt(abs(sum((x-x0).^2 .* abs(y2)) / abs(sum(y2))));
        sigma_est = 0.3;
        %sigma_est = (max(x) - min(x)) / 20;
        x0 = [
            A0
            A1
            A2
            sigma_est
            sigma_est
            sigma_est
            Delta
            x0
            bg
        ];


    end

    %% options
    % the base_fit function is designed to be flexible, and can accept a
    % lot of different parameters to adjust 
    %%%%% Default Options %%%%%
%     options = struct(...
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
    
% coeffs = [A1, A0, A2, sigma1, sigma0, sigma2, Delta, x0, bg]
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', true,...
        'PlotInitialGuess', true,...
        'XAxisLabel', 'mmWave WindFreak Synth (MHz)' ,...
        'YAxisLabel', 'MCS Signal Ratio',...
        'CoeffNames', {{'Area 1', 'Area 2', 'Area 3', 'FWHM 1', 'FWHM 2', 'FWHM 3', 'Zeeman Splitting', 'Center Peak','Offset'}},...
        'CoeffUnits', {{'','','','MHz','MHz','MHz','MHz','MHz','',''}},...
        'AnnotateFunction', @myAnnotate,...
        'FitTitle', 'Triple Lorentzian Fit',...
        'Statistics', 'gaussian');
    
    [xav,yav,yer,coefflist,coefflist_err] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @initial_guess, options)
    
    % if the saveVals is toggled on, save the values from the data
    
    if (saveVals == 1)
        filename = "MMWavePlotVals.xlsx";
        sheet    = "Sheet1";

        timeStr = strjoin(string(analyVar.timevectorAtom), ",");   % "1.2,3.4,5.6"
        plugInStr = strjoin(string(analyVar.plugInVec), ",");
        
        % Save several values to the array
        %first calculate the polarization purity in function too... (ctrl
        %+c ctrl + v from the myAnnotate function)
        % evaluating the integral from the fit coefficients.
        lineintegral1 = coefflist{1}(1); % amp*FWHM*Pi/2 = lorentzian integral
        lineintegral2 = coefflist{1}(2); % amp*FWHM*Pi/2 = lorentzian integral
        lineintegral3 = coefflist{1}(3); % amp*FWHM*Pi/2 = lorentzian integral
        %lineintegral_err1 = lineintegral1*(coefflist_err{1}(1)/coefflist{1}(1)+ coefflist{1}(4)/coefflist{1}(4));
        %lineintegral_err2 = lineintegral2*(coefflist_err{1}(2)/coefflist{1}(2)+ coefflist{1}(5)/coefflist{1}(5));
        %lineintegral_err3 = lineintegral3*(coefflist{1}(3)/coefflist{1}(3)+ coefflist{1}(6)/coefflist{1}(6));
    
        %Ratio of purity
        purityRatio = abs((lineintegral2)/(lineintegral2+lineintegral1+lineintegral3));
        purityError = (sqrt((coefflist_err{1}(2)/coefflist{1}(2))^2+(coefflist_err{1}(1)+coefflist_err{1}(2) + coefflist_err{1}(3)/coefflist{1}(1) + coefflist{1}(2) + coefflist{1}(3))^2)) ;
        %%%
        % time stamp | plug in vector | Amplitude 1 | Amplitude 2 |
        % Amplitude 3 | FWHM 1 | FWHM 2 | FWHM 3 | Zeeman Splitting |
        % Center Peak | Offset | Purity Ratio | first file name
        row = {timeStr, plugInStr, coefflist{1}(1), coefflist{1}(2), coefflist{1}(3), coefflist{1}(4), coefflist{1}(5), coefflist{1}(6), coefflist{1}(7), coefflist{1}(8), coefflist{1}(9), purityRatio,'58 s', '58 d', indivDataset{1,1}.fileAtom{1}};  % cell row
        
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
function h = myplot(x,y,analyVar,i)
    h = plot(x,y,...
        'LineStyle','none',...
        'Marker', 'o',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(i,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(i,:));
end

function h = myerrorbar(x,y,yerr,analyVar,i)
    h = errorbar(x,y,yerr,...
        'LineStyle','none',...
        'Marker', 'o',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(i,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(i,:));
end

function h = myfitplot(x,y,analyVar,i)
    h = plot(x,y,...
        'LineStyle','-',...
        'Marker', 'none',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(i,:),...
        'MarkerEdgeColor', 'none',...
        'Color', analyVar.COLORS(i,:),...
        'LineWidth',5);
end

function an = myAnnotate(coeffs, err, coeffNames, coeffUnits)

    % modify coefficients to be in the correct units

    % make HWHM FWHM
    coeffs(4) = coeffs(4)*2;
    coeffs(5) = coeffs(5)*2;
    coeffs(6) = coeffs(6)*2;

    dim = [.6 .5 .3 .3];

    %% Coeff(1-3) are fitted to be the Area of the Lorentzian.
    %% Area_i ~ E_i where i is the peak and E is the polarization vector
    %% Then Pi Purity Ratio should be E_i / Sqrt[ Sum(E_i^2) ], the projection of the E-field on the z-axis

    %Ratio of purity
    % purityRatio = abs((coeffs(2))/(coeffs(2) + coeffs(3) + coeffs(1)));
    % purityError = sqrt((err(2)/coeffs(2))^2+((err(1)+err(2) + err(3))/(coeffs(1) + coeffs(2) + coeffs(3)))^2) * purityRatio;
    purityRatio = abs((coeffs(2))/sqrt(coeffs(2)^2 + coeffs(3)^2 + coeffs(1)^2));
    purityError = sqrt( ((coeffs(1)^2+coeffs(3)^2)^2*err(2)^2 + (coeffs(2)*coeffs(1)*err(1))^2 + (coeffs(2)*coeffs(3)*err(3))^2 ) / ...
        (coeffs(1)^2 + coeffs(2)^2 + coeffs(3)^2)^3 ) * purityRatio;
    %%%
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
    temp = .001;
    strs = cell(numel(coeffs)+1,1);
    for i = 1:numel(coeffs)
        if i < numel(coeffs)
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}, newline];
        else
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}];
        end
    end
    %%% adding line integral string.
    strs{i+1} = [newline,'PurityPiRatio', ': ', num2str(purityRatio)];
    strs{i+2} = [newline,'Purity Error', ': ', num2str(purityError)];

    
    an = annotation('textbox', dim, 'String', strjoin(strs),...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
end