function funcOut = zeeman_sfi_triple_gaussian_loss(analyVar, indivDataset, avgDataset)
    
    % Function that fits to a triple-peaked Gaussian signal-loss spectrum.
    %
    % The three features are modeled as TROUGHS on a higher background:
    %
    % y = bg ...
    %     - A1*exp(-(x-(x0-Delta))^2/(2*sigma1^2)) ...
    %     - A0*exp(-(x-x0)^2/(2*sigma0^2)) ...
    %     - A2*exp(-(x-(x0+Delta))^2/(2*sigma2^2))
    %
    % coeffs = [A1, A0, A2, sigma1, sigma0, sigma2, Delta, x0, bg]
    %
    % A1/A0/A2 are positive signal-loss amplitudes.
    % sigma1/sigma0/sigma2 are Gaussian sigma widths.
    
    form = @(coeffs,x) ...
        coeffs(9) ...
        - coeffs(1)*exp(-(x-(coeffs(8)-coeffs(7))).^2 ./ (2*coeffs(4)^2)) ...
        - coeffs(2)*exp(-(x-coeffs(8)).^2 ./ (2*coeffs(5)^2)) ...
        - coeffs(3)*exp(-(x-(coeffs(8)+coeffs(7))).^2 ./ (2*coeffs(6)^2));

    indVarField = 'imagevcoAtom';
    depVarField = 'OD_TotalCounts';
    
    %% toggle saveVals on and off to save values to an external excel sheet
    saveVals = 0;
   
    function x0 = initial_guess(x, y)

        % Sort x values in ascending order
        [xs, idx] = sort(x);
        ys = y(idx);

        % ---------------------------------------------------------
        % Background
        % ---------------------------------------------------------
        bg = median(ys(1:round(end*0.1))); 
        bg = 0;
        
        % Signal loss = background - measured signal
        % This should be positive at the troughs.
        yLoss = bg + ys;
        %yLoss = bg;

        % ---------------------------------------------------------
        % Estimate center peak position
        % ---------------------------------------------------------
        % Weight by signal loss rather than signal itself
        if sum(yLoss) > 0
            centerGuess = sum(xs .* yLoss) / sum(yLoss);
        else
            centerGuess = mean(xs);
        end

        centerGuess = 81.93;

        % x0 = 12897;
        % x0 = 235.4;

        % ---------------------------------------------------------
        % Estimate Zeeman splitting
        % ---------------------------------------------------------
        Delta = 0.08;
        disp(Delta)

        % ---------------------------------------------------------
        % Amplitudes = depth of signal-loss troughs
        % ---------------------------------------------------------
        A1 = interp1(xs, yLoss, centerGuess-Delta, ...
            'linear', 'extrap');

        A0 = interp1(xs, yLoss, centerGuess, ...
            'linear', 'extrap');

        A2 = interp1(xs, yLoss, centerGuess+Delta, ...
            'linear', 'extrap');

        % Make sure amplitudes are positive
        A1 = max(A1, 0.001);
        A0 = max(A0, 0.001);
        A2 = max(A2, 0.001);
        A1 = -800;
        A0 = -800;
        A2 = -800;

        % ---------------------------------------------------------
        % Widths
        % ---------------------------------------------------------
        sigma_est = 0.04;

        % ---------------------------------------------------------
        % Initial coefficient vector
        %
        % [A1, A0, A2, sigma1, sigma0, sigma2, Delta, x0, bg]
        % ---------------------------------------------------------
        x0 = [
            A1
            A0
            A2
            sigma_est
            sigma_est
            sigma_est
            Delta
            centerGuess
            bg
        ];

    end

    %% options
    
    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', true,...
        'PlotInitialGuess', true,...
        'XAxisLabel', '689nm rMOT Frequency (MHz)' ,...
        'YAxisLabel', 'Integrated Background Subtracted Counts',...
        'CoeffNames', {{ ...
            'Area 1', ...
            'Area 2', ...
            'Area 3', ...
            'FWHM 1', ...
            'FWHM 2', ...
            'FWHM 3', ...
            'Zeeman Splitting', ...
            'Center Peak', ...
            'Offset'}},...
        'CoeffUnits', {{ ...
            '', ...
            '', ...
            '', ...
            'MHz', ...
            'MHz', ...
            'MHz', ...
            'MHz', ...
            'MHz', ...
            ''}},...
        'AnnotateFunction', @myAnnotate,...
        'FitTitle', 'Triple Gaussian Signal Loss Fit',...
        'Statistics', 'gaussian');
    
    [xav,yav,yer,coefflist,coefflist_err] = ...
        base_fit(analyVar, indivDataset, avgDataset, ...
        form, indVarField, depVarField, @initial_guess, options);
    
    %% Save values to external Excel sheet
    
    if (saveVals == 1)

        filename = "MMWavePlotVals.xlsx";
        sheet = "Sheet1";

        timeStr = strjoin(string(analyVar.timevectorAtom), ",");
        plugInStr = strjoin(string(analyVar.plugInVec), ",");

        % ---------------------------------------------------------
        % Gaussian integrated signal-loss areas
        %
        % Integral = A*sigma*sqrt(2*pi)
        % ---------------------------------------------------------
        
        lineintegral1 = coefflist{1}(1) * ...
            coefflist{1}(4) * sqrt(2*pi);

        lineintegral2 = coefflist{1}(2) * ...
            coefflist{1}(5) * sqrt(2*pi);

        lineintegral3 = coefflist{1}(3) * ...
            coefflist{1}(6) * sqrt(2*pi);

        % ---------------------------------------------------------
        % Purity ratio
        %
        % Since Gaussian areas depend on both amplitude and width,
        % use the integrated areas rather than amplitudes alone.
        % ---------------------------------------------------------
        
        purityRatio = abs(lineintegral2 / ...
            sqrt(lineintegral2^2 + ...
                 lineintegral1^2 + ...
                 lineintegral3^2));

        % ---------------------------------------------------------
        % Uncertainty in integrated areas
        % ---------------------------------------------------------
        
        lineintegral1_err = abs(lineintegral1) * sqrt( ...
            (coefflist_err{1}(1)/coefflist{1}(1))^2 + ...
            (coefflist_err{1}(4)/coefflist{1}(4))^2);

        lineintegral2_err = abs(lineintegral2) * sqrt( ...
            (coefflist_err{1}(2)/coefflist{1}(2))^2 + ...
            (coefflist_err{1}(5)/coefflist{1}(5))^2);

        lineintegral3_err = abs(lineintegral3) * sqrt( ...
            (coefflist_err{1}(3)/coefflist{1}(3))^2 + ...
            (coefflist_err{1}(6)/coefflist{1}(6))^2);

        % ---------------------------------------------------------
        % Purity uncertainty
        % ---------------------------------------------------------
        
        denom = sqrt(lineintegral1^2 + ...
                     lineintegral2^2 + ...
                     lineintegral3^2);

        purityError = sqrt( ...
            ((lineintegral1^2 + lineintegral3^2)^2 * ...
                lineintegral2_err^2 + ...
            (lineintegral2*lineintegral1*lineintegral1_err)^2 + ...
            (lineintegral2*lineintegral3*lineintegral3_err)^2) ...
            / denom^6 );

        % ---------------------------------------------------------
        % Convert sigma -> FWHM for saved values
        % ---------------------------------------------------------
        
        fwhm1 = coefflist{1}(4) * 2*sqrt(2*log(2));
        fwhm2 = coefflist{1}(5) * 2*sqrt(2*log(2));
        fwhm3 = coefflist{1}(6) * 2*sqrt(2*log(2));

        % ---------------------------------------------------------
        % Row:
        %
        % time stamp | plug in vector | Loss 1 | Loss 2 |
        % Loss 3 | FWHM 1 | FWHM 2 | FWHM 3 |
        % Zeeman Splitting | Center Peak | Offset |
        % Purity Ratio | first file name
        % ---------------------------------------------------------
        
        row = { ...
            timeStr, ...
            plugInStr, ...
            coefflist{1}(1), ...
            coefflist{1}(2), ...
            coefflist{1}(3), ...
            fwhm1, ...
            fwhm2, ...
            fwhm3, ...
            coefflist{1}(7), ...
            coefflist{1}(8), ...
            coefflist{1}(9), ...
            purityRatio, ...
            '58 s', ...
            '58 d', ...
            indivDataset{1,1}.fileAtom{1}};

        if isfile(filename)
            C = readcell(filename, "Sheet", sheet);
            nextRow = size(C,1) + 1;
        else
            nextRow = 1;
        end

        range = "A" + nextRow;

        writecell(row, filename, ...
            "Sheet", sheet, ...
            "Range", range);

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

    dim = [.6 .5 .3 .3];

    % ---------------------------------------------------------
    % Convert Gaussian sigma -> FWHM
    % ---------------------------------------------------------
    
    coeffs(4) = coeffs(4) * 2*sqrt(2*log(2));
    coeffs(5) = coeffs(5) * 2*sqrt(2*log(2));
    coeffs(6) = coeffs(6) * 2*sqrt(2*log(2));

    err(4) = err(4) * 2*sqrt(2*log(2));
    err(5) = err(5) * 2*sqrt(2*log(2));
    err(6) = err(6) * 2*sqrt(2*log(2));

    % ---------------------------------------------------------
    % Gaussian integrated signal-loss areas
    %
    % A*sigma*sqrt(2*pi)
    %
    % Need sigma, so convert FWHM back to sigma.
    % ---------------------------------------------------------
    
    sigma1 = coeffs(4) / (2*sqrt(2*log(2)));
    sigma2 = coeffs(5) / (2*sqrt(2*log(2)));
    sigma3 = coeffs(6) / (2*sqrt(2*log(2)));

    sigma1_err = err(4) / (2*sqrt(2*log(2)));
    sigma2_err = err(5) / (2*sqrt(2*log(2)));
    sigma3_err = err(6) / (2*sqrt(2*log(2)));

    area1 = coeffs(1) * sigma1 * sqrt(2*pi);
    area2 = coeffs(2) * sigma2 * sqrt(2*pi);
    area3 = coeffs(3) * sigma3 * sqrt(2*pi);

    % ---------------------------------------------------------
    % Area uncertainties
    % ---------------------------------------------------------
    
    area1_err = abs(area1) * sqrt( ...
        (err(1)/coeffs(1))^2 + ...
        (sigma1_err/sigma1)^2);

    area2_err = abs(area2) * sqrt( ...
        (err(2)/coeffs(2))^2 + ...
        (sigma2_err/sigma2)^2);

    area3_err = abs(area3) * sqrt( ...
        (err(3)/coeffs(3))^2 + ...
        (sigma3_err/sigma3)^2);

    % ---------------------------------------------------------
    % Purity ratio
    %
    % Use integrated Gaussian areas because the area depends
    % on both amplitude and width.
    % ---------------------------------------------------------
    
    denom = sqrt(area1^2 + area2^2 + area3^2);

    purityRatio = abs(area2 / denom);

    % ---------------------------------------------------------
    % Propagate purity uncertainty
    % ---------------------------------------------------------
    
    purityError = sqrt( ...
        ((area1^2 + area3^2)^2 * area2_err^2 + ...
        (area2*area1*area1_err)^2 + ...
        (area2*area3*area3_err)^2) / ...
        denom^6);

    % ---------------------------------------------------------
    % Default names/units
    % ---------------------------------------------------------
    
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

    % ---------------------------------------------------------
    % Build annotation
    % ---------------------------------------------------------
    
    strs = cell(numel(coeffs)+3,1);

    for i = 1:numel(coeffs)

        if i < numel(coeffs)

            strs{i} = [ ...
                coeffNames{i}, ': ', ...
                unc_string(coeffs(i),err(i)), ...
                ' ', coeffUnits{i}, newline];

        else

            strs{i} = [ ...
                coeffNames{i}, ': ', ...
                unc_string(coeffs(i),err(i)), ...
                ' ', coeffUnits{i}];

        end

    end

    % Add purity information
    strs{numel(coeffs)+1} = [ ...
        newline, ...
        'PurityPiRatio: ', ...
        num2str(purityRatio)];

    strs{numel(coeffs)+2} = [ ...
        newline, ...
        'Purity Error: ', ...
        num2str(purityError)];

    strs{numel(coeffs)+3} = [ ...
        newline, ...
        'Center Freq Error: ', ...
        num2str(err(8))];

    an = annotation('textbox', dim, ...
        'String', strjoin(strs),...
        'FitBoxToText', 'on',...
        'BackgroundColor', 'white');

end