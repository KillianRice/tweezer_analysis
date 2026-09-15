function funcOut = skewed_sfi_gaussian(analyVar, indivDataset, avgDataset)

    % =========================================================
    % FIT MODE
    %
    % fitMode = 1 --> Single skewed Gaussian
    % fitMode = 2 --> Double skewed Gaussian
    % fitMode = 3 --> Normal Single Gaussian
    % =========================================================

    fitMode = 1;

    % Typical parameters for double skewed Gaussian
    typicalFWHM = 0.2;       % MHz
    typicalSpacing = 0.1;    % MHz

    % ---------------------------------------------------------
    % Single Gaussian
    % ---------------------------------------------------------

    form = @(coeffs,x) ...
        coeffs(1) * exp(-(x-coeffs(2)).^2 ./ ...
        (2*coeffs(3)^2)) + coeffs(4);

    % ---------------------------------------------------------
    % Single skewed Gaussian
    % coeffs = [A, mu, sigma, offset, alpha]
    % ---------------------------------------------------------

    formSkewed = @(coeffs,x) ...
        (coeffs(1) * ...
        exp(-(x-coeffs(2)).^2 ./ (2*coeffs(3)^2)) .* ...
        2 .* normcdf(coeffs(5) * ...
        (x-coeffs(2)) / coeffs(3))) + coeffs(4);

    % ---------------------------------------------------------
    % Double skewed Gaussian
    %
    % coeffs =
    %
    % [A1, A2, sigma1, sigma2, Delta, x0, alpha1, alpha2, offset]
    %
    % Peak 1 center = x0
    % Peak 2 center = x0 + Delta
    %
    % For signal loss:
    % A1 and A2 should be NEGATIVE.
    % ---------------------------------------------------------

    formDoubleSkewed = @(coeffs,x) ...
        coeffs(9) + ...
        coeffs(1) .* ...
        exp(-(x-coeffs(6)).^2 ./ ...
        (2*coeffs(3)^2)) .* ...
        2 .* normcdf(coeffs(7) .* ...
        (x-coeffs(6)) ./ coeffs(3)) + ...
        coeffs(2) .* ...
        exp(-(x-(coeffs(6)+coeffs(5))).^2 ./ ...
        (2*coeffs(4)^2)) .* ...
        2 .* normcdf(coeffs(8) .* ...
        (x-(coeffs(6)+coeffs(5))) ./ coeffs(4));

    indVarField = 'imagevcoAtom';

    % depVarField = 'sfiIntegral_roi1';
    %depVarField = 'sfiIntegral';

    depVarField = 'OD_TotalCounts';


    % =========================================================
    % INITIAL GUESS: SINGLE SKEWED GAUSSIAN
    % =========================================================

    function x0 = initial_guess_skewed(x, y)

        x0 = zeros(5,1);

        x = x(:);
        y = y(:);

        % Baseline
        x0(4) = min(y);

        % Negative amplitude for signal loss
        x0(1) = max(y);

        % Inverted profile
        y_inv = x0(4) - y;

        total_area = sum(y_inv);

        if total_area == 0
            total_area = 1;
        end

        % Center
        com_guess = sum(x .* y_inv) / total_area;

        % Sigma
        sigma_guess = sqrt( ...
            sum((x-com_guess).^2 .* y_inv) / total_area);

        if sigma_guess == 0
            sigma_guess = 1;
        end

        % Your experimentally determined sigma
        sigma_guess = 0.1;

        % Initial skew
        x0(5) = 2;

        % Correct mu for skewness
        delta_guess = ...
            x0(5) / sqrt(1+x0(5)^2);

        x0(2) = com_guess - ...
            sigma_guess * delta_guess * sqrt(2/pi);

        x0(3) = sigma_guess;

    end


    % =========================================================
    % INITIAL GUESS: NORMAL GAUSSIAN
    % =========================================================

    function x0 = initial_guess(x,y)

        x0 = zeros(4,1);
        %x0(1) = (min(y)-max(y)); %% For negative gaussians
        x0(1) = max(y);             %% For positive
        x0(2) =  sum(x.*y)/sum(y);
        %x0(2) = 81.96;
        x0(3) = sqrt(sum((x-x0(2)).^2.*y)/sum(y));
        x0(3) = .2;
        %x0(4) = max(y);            %% For negative gaussians
        x0(4) = min(y);             %% For positive

    end


% =========================================================
% INITIAL GUESS: DOUBLE SKEWED GAUSSIAN
%
% peakSign = -1  -> negative peaks (signal loss)
% peakSign = +1  -> positive peaks
% =========================================================

function x0 = initial_guess_double_skewed(x,y)

    x = x(:);
    y = y(:);

    % -----------------------------------------------------
    % Peak sign
    % -----------------------------------------------------

    peakSign = 1;   % Default: signal loss

    if ~ismember(peakSign,[-1 1])
        error('peakSign must be +1 or -1.');
    end

    % -----------------------------------------------------
    % Background
    % -----------------------------------------------------

    nEdge = max(3,round(0.1*numel(y)));

    offset = median([ ...
        y(1:nEdge); ...
        y(end-nEdge+1:end)]);

    % -----------------------------------------------------
    % Convert signal into positive peak strength
    %
    % peakSign = -1:
    %   negative peaks -> positive loss
    %
    % peakSign = +1:
    %   positive peaks -> positive signal
    % -----------------------------------------------------

    yPeak = peakSign .* (y - offset);

    yPeak = max(yPeak,0);

    % -----------------------------------------------------
    % Find two strongest peaks
    % -----------------------------------------------------

    try

        [pks,locs,widths] = findpeaks( ...
            yPeak,...
            x,...
            'MinPeakDistance',0.05,...
            'SortStr','descend');

    catch

        pks = [];
        locs = [];
        widths = [];

    end

    % -----------------------------------------------------
    % TWO PEAKS FOUND
    % -----------------------------------------------------

    if numel(locs) >= 2

        % Keep two strongest
        pks  = pks(1:2);
        locs = locs(1:2);

        % Sort low frequency -> high frequency
        [locs,order] = sort(locs);

        pks = pks(order);

        % -------------------------------------------------
        % Amplitudes
        % -------------------------------------------------

        A1 = peakSign * max(pks(1),0.001);
        A2 = peakSign * max(pks(2),0.001);

        % -------------------------------------------------
        % Separation
        % -------------------------------------------------

        Delta = locs(2) - locs(1);

        % -------------------------------------------------
        % Center of first peak
        % -------------------------------------------------

        xCenter = locs(1);

        % -------------------------------------------------
        % Widths
        % -------------------------------------------------

        if numel(widths) >= 2

            sigma1 = widths(1) / ...
                (2*sqrt(2*log(2)));

            sigma2 = widths(2) / ...
                (2*sqrt(2*log(2)));

        else

            sigma1 = typicalFWHM / ...
                (2*sqrt(2*log(2)));

            sigma2 = typicalFWHM / ...
                (2*sqrt(2*log(2)));

        end

    else

        % =================================================
        % FALLBACK
        % =================================================

        disp('Could not identify two peaks.');
        disp('Using typical double-skewed Gaussian guess.');

        [maxPeak,maxIdx] = max(yPeak);

        xCenter = x(maxIdx);

        A1 = peakSign * max(maxPeak,0.001);
        A2 = peakSign * max(maxPeak,0.001);

        Delta = typicalSpacing;

        sigma1 = typicalFWHM / ...
            (2*sqrt(2*log(2)));

        sigma2 = typicalFWHM / ...
            (2*sqrt(2*log(2)));

    end

    % -----------------------------------------------------
    % Initial skew parameters
    % -----------------------------------------------------

    alpha1 = 2;
    alpha2 = 2;

    % sigma2 = sigma2 * 10;
    % A1 = 600;
    % A2 = 200;
    % Delta = 0.6;

    % -----------------------------------------------------
    % Coefficient vector
    %
    % [A1,A2,sigma1,sigma2,Delta,x0,alpha1,alpha2,offset]
    % -----------------------------------------------------

    x0 = [ ...
        A1
        A2
        sigma1
        sigma2
        Delta
        xCenter
        alpha1
        alpha2
        offset];

end



    %% ========================================================
    % OPTIONS
    % =========================================================

    if fitMode == 1

        % -----------------------------------------------------
        % SINGLE SKEWED GAUSSIAN
        % -----------------------------------------------------

        options = struct(...
            'PlotIndivFits',false,...
            'PlotAll',true,...
            'PlotAllAvgs',true,...
            'PlotInitialGuess',true,...
            'XAxisLabel','Vertical Red Cooling Frequency (MHz)',...
            'CoeffNames',{{ ...
                'Amplitude',...
                'Line Center',...
                'FWHM[For Symm. Gaus.]',...
                'Offset',...
                'alpha'}},...
            'CoeffUnits',{{ ...
                '',...
                'MHz',...
                'MHz',...
                '',...
                ''}},...
            'AnnotateFunction',@myAnnotateSkewed,...
            'Statistics','gaussian',...
            'PlotNormalized',false,...
            'PlotNormToMeasuredParams',false);

        [xav,yav,yer,coefflist,coefflist_err,...
            avg_fit_coeffs_twzr] = ...
            base_fit(analyVar,...
            indivDataset,...
            avgDataset,...
            formSkewed,...
            indVarField,...
            depVarField,...
            @initial_guess_skewed,...
            options);


    elseif fitMode == 2

        % -----------------------------------------------------
        % DOUBLE SKEWED GAUSSIAN
        % -----------------------------------------------------

        options = struct(...
            'PlotIndivFits',false,...
            'PlotAll',true,...
            'PlotAllAvgs',true,...
            'PlotInitialGuess',true,...
            'XAxisLabel','Valon Synth (MHz)',...
            'CoeffNames',{{ ...
                'Amplitude 1',...
                'Amplitude 2',...
                'FWHM 1',...
                'FWHM 2',...
                'Peak Separation',...
                'Line Center 1',...
                'alpha 1',...
                'alpha 2',...
                'Offset'}},...
            'CoeffUnits',{{ ...
                '',...
                '',...
                'MHz',...
                'MHz',...
                'MHz',...
                'MHz',...
                '',...
                '',...
                ''}},...
            'AnnotateFunction',@myAnnotateDoubleSkewed,...
            'Statistics','gaussian',...
            'PlotNormalized',false,...
            'PlotNormToMeasuredParams',false);

        [xav,yav,yer,coefflist,coefflist_err,...
            avg_fit_coeffs_twzr] = ...
            base_fit(analyVar,...
            indivDataset,...
            avgDataset,...
            formDoubleSkewed,...
            indVarField,...
            depVarField,...
            @initial_guess_double_skewed,...
            options);

    elseif fitMode == 3

        % -----------------------------------------------------
        % SINGLE UNSKEWED GAUSSIAN
        % -----------------------------------------------------

        options = struct(...
            'PlotIndivFits',false,...
            'PlotAll',true,...
            'PlotAllAvgs',true,...
            'PlotInitialGuess',true,...
            'XAxisLabel','Valon Synth (MHz)',...
            'CoeffNames',{{ ...
                'Amplitude',...
                'Line Center',...
                'FWHM',...
                'Offset'}},...
            'CoeffUnits',{{ ...
                '',...
                'MHz',...
                'MHz',...
                ''}},...
            'AnnotateFunction',@myAnnotate,...
            'Statistics','gaussian',...
            'PlotNormalized',false,...
            'PlotNormToMeasuredParams',false);

        [xav,yav,yer,coefflist,coefflist_err,...
            avg_fit_coeffs_twzr] = ...
            base_fit(analyVar,...
            indivDataset,...
            avgDataset,...
            form,...
            indVarField,...
            depVarField,...
            @initial_guess,...
            options);

    end


    avgDataset.avg_fit_coeffs_twzr = avg_fit_coeffs_twzr;


    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end

% function funcOut = sfi_gaussian(analyVar, indivDataset, avgDataset)
% 
%     form = @(coeffs,x) coeffs(1) * exp(-(x-coeffs(2)).^2 ./ (2*coeffs(3)^2)) + coeffs(4);
% 
%     formSkewed = @(coeffs,x) (coeffs(1) * exp(-(x-coeffs(2)).^2 ./ (2*coeffs(3)^2)) .* 2 .* normcdf(coeffs(5) * (x - coeffs(2)) / coeffs(3))) + coeffs(4); 
% 
%     fitSkew = 1; %% Option for fitting to a skewed Gaussian... 
% 
%     indVarField = 'imagevcoAtom';
%     % depVarField = 'sfiIntegral_roi1';
%     %depVarField = 'sfiIntegral';
%     depVarField = 'OD_TotalCounts';
% 
%     function x0 = initial_guess_skewed(x, y)
%         x0 = zeros(5,1);
% 
%         % Ensure x and y are column vectors for vector math
%         x = x(:);
%         y = y(:);
% 
%         % 1. Baseline offset (Top level for an inverted peak)
%         x0(4) = max(y); 
% 
%         % 2. Amplitude (Negative value because it dips down from baseline)
%         x0(1) = min(y) - x0(4); 
% 
%         % Create an inverted, positive-facing profile for moment calculations
%         y_inv = x0(4) - y;
%         total_area = sum(y_inv);
% 
%         % Prevent divide-by-zero if data is completely flat
%         if total_area == 0
%             total_area = 1;
%         end
% 
%         % 3. True Center of Mass (Mean of the inverted peak)
%         com_guess = sum(x .* y_inv) / total_area;
%         com_guess = 82;
% 
%         % 4. Empirical Sigma (Standard deviation of the inverted peak)
%         sigma_guess = sqrt(sum((x - com_guess).^2 .* y_inv) / total_area);
%         if sigma_guess == 0, sigma_guess = 1; end % Fallback safety
% 
%         sigma_guess = 0.1;
% 
%         % 5. Initial Skewness Guess
%         % A value of 2 provides a moderate right-hand tail bias to start the fit
%         x0(5) = 2; 
% 
%         % 6. Location Parameter (mu)
%         % Since skewness is positive (>0), mu must sit to the LEFT of the center of mass.
%         delta_guess = x0(5) / sqrt(1 + x0(5)^2);
%         x0(2) = com_guess - sigma_guess * delta_guess * sqrt(2/pi);
% 
%         % 7. Sigma parameter assignment
%         x0(3) = sigma_guess;
%     end
% 
% 
%     function x0 = initial_guess(x,y)
%             x0 = zeros(4,1);
%             x0(1) = (min(y)-max(y)); %% For negative gaussians
%             %x0(1) = max(y);             %% For positive
%             %x0(2) =  sum(x.*y)/sum(y);
%             x0(2) = 82.1;
%             x0(3) = sqrt(sum((x-x0(2)).^2.*y)/sum(y));
%             x0(3) = .02;
%             x0(4) = max(y);            %% For negative gaussians
%             %x0(4) = min(y);             %% For positive
%     end
% 
%     %% options
%     % the base_fit function is designed to be flexible, and can accept a
%     % lot of different parameters to adjust 
%     %%%%% Default Options %%%%%
% %     options = struct(...
% %         'DataPlotFunction', @defaultDataPlot,...
% %         'AvgDataPlotFunction', @defaultAvgDataPlot,...
% %         'FitLinePlotFunction', @defaultFitLinePlot,...
% %         'AnnotateFunction', @defaultAnnotate,...
% %         'IndivFitPlotFunction', @defaultIndivFitPlot,...
% %         'PlotIndivFits' , true,...
% %         'PlotAvgFits', true,...
% %         'XAxisLabel', indVarField ,...
% %         'YAxisLabel', depVarField,...
% %         'FitLB', [],...
% %         'FitUB', [],...
% %         'FitOptions', struct('Display','off'),...
% %         'PlotInitialGuess', true, ...
% %         'InitialGuessPlotFunction', @defaultInitialGuessPlot,...
% %         'PlotAll', true,...
% %         'PlotAllAvgs', true, ...
% %         'CoeffNames', {{}},...
% %         'CoeffUnits', {{}},...
% %         'YAxisScale', 'linear',...
% %         'XAxisScale', 'linear');
% 
% 
%     if fitSkew ==1 
%         options = struct(...
%         'PlotIndivFits', false,...
%         'PlotAll', true,...
%         'PlotAllAvgs', true,...
%         'PlotInitialGuess', true,...
%         'XAxisLabel', 'Valon Synth (MHz)',...
%         'CoeffNames', {{'Amplitude', 'Line Center', 'FWHM[For Symm. Gaus.]', 'Offset', 'alpha'}},...
%         'CoeffUnits', {{'','MHz','MHz','',''}},...
%         'AnnotateFunction', @myAnnotateSkewed,...
%         'Statistics', 'gaussian', ...
%         'PlotNormalized', false, ...
%         'PlotNormToMeasuredParams', false);
% 
%         [xav,yav,yer,coefflist,coefflist_err,avg_fit_coeffs_twzr] = base_fit(analyVar, indivDataset, avgDataset, formSkewed, indVarField, depVarField, @initial_guess_skewed, options);
% 
% 
%     else
%         options = struct(...
%         'PlotIndivFits', false,...
%         'PlotAll', true,...
%         'PlotAllAvgs', true,...
%         'PlotInitialGuess', false,...
%         'XAxisLabel', 'Valon Synth (MHz)',...
%         'CoeffNames', {{'Amplitude', 'Line Center', 'FWHM', 'Offset'}},...
%         'CoeffUnits', {{'','MHz','MHz',''}},...
%         'AnnotateFunction', @myAnnotate,...
%         'Statistics', 'gaussian', ...
%         'PlotNormalized', true, ...
%         'PlotNormToMeasuredParams', false);
%         [~,~,~,~,~,avg_fit_coeffs_twzr] = base_fit(analyVar, indivDataset, avgDataset, form, indVarField, depVarField, @initial_guess, options);
%     end
% 
%     avgDataset.avg_fit_coeffs_twzr = avg_fit_coeffs_twzr;
% 
% 
%     funcOut.analyVar = analyVar;
%     funcOut.indivDataset = indivDataset;
%     funcOut.avgDataset = avgDataset;
% 
% end
% function h = myplot(x,y,analyVar,i)
%     h = plot(x,y,...
%         'LineStyle','none',...
%         'Marker', 'o',...
%         'MarkerSize', analyVar.markerSize,...
%         'MarkerFaceColor', analyVar.COLORS(i,:),...
%         'MarkerEdgeColor', 'none',...
%         'Color', analyVar.COLORS(i,:));
% end
% 
% function h = myerrorbar(x,y,yerr,analyVar,i)
%     h = errorbar(x,y,yerr,...
%         'LineStyle','none',...
%         'Marker', 'o',...
%         'MarkerSize', analyVar.markerSize,...
%         'MarkerFaceColor', analyVar.COLORS(i,:),...
%         'MarkerEdgeColor', 'none',...
%         'Color', analyVar.COLORS(i,:));
% end
% 
% function h = myfitplot(x,y,analyVar,i)
%     h = plot(x,y,...
%         'LineStyle','-',...
%         'Marker', 'none',...
%         'MarkerSize', analyVar.markerSize,...
%         'MarkerFaceColor', analyVar.COLORS(i,:),...
%         'MarkerEdgeColor', 'none',...
%         'Color', analyVar.COLORS(i,:),...
%         'LineWidth',5);
% end
% 
function an = myAnnotate(coeffs, err, coeffNames, coeffUnits)

    dim = [.6 .5 .3 .3];

    % convert to relevant parameters
    coeffs(3) = coeffs(3) * 2*sqrt(2*log(2));
    err(3) = err(3) * 2*sqrt(2*log(2));
    % evaluating the integral from the fit coefficients.
    lineintegral = coeffs(1)*coeffs(3)*sqrt(2*pi); % amp*FWHM*Sqrt(2*Pi) = Gaussian integral
    lineintegral_err = lineintegral*(err(3)/coeffs(3)+ err(1)/coeffs(1));
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

    strs = cell(numel(coeffs)+1,1);
    for i = 1:numel(coeffs)
        if i < numel(coeffs)
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i},newline];
        else
            strs{i} = [coeffNames{i}, ': ', unc_string(coeffs(i),err(i)),...
                ' ', coeffUnits{i}];
        end
    end
    %%% adding line integral string.
    strs{i+1} = [newline,'LineIntegral', ': ', unc_string(lineintegral,lineintegral_err)];

    an = annotation('textbox', dim, 'String', strjoin(strs),...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
end

function an = myAnnotateSkewed(coeffs, err, coeffNames, coeffUnits)

    dim = [.6 .5 .3 .3];

    % ---------------------------------------------------------
    % Skewed Gaussian parameters
    %
    % coeffs(1) = amplitude
    % coeffs(2) = location parameter (mu)
    % coeffs(3) = sigma
    % coeffs(4) = offset
    % coeffs(5) = alpha
    % ---------------------------------------------------------

    sigma = coeffs(3);
    sigma_err = err(3);

    % Convert sigma -> FWHM for display
    fwhm_sym = sigma * 2*sqrt(2*log(2));
    fwhm_sym_err = sigma_err * 2*sqrt(2*log(2));

    % ---------------------------------------------------------
    % Reconstruct Fitted Function & Solve for Empirical FWHM
    % ---------------------------------------------------------
    % 1. Rebuild the anonymous function handle using the local coeffs array
    formSkewed_local = @(c, x) (c(1) * exp(-(x - c(2)).^2 ./ (2 * c(3)^2)) .* 2 .* normcdf(c(5) * (x - c(2)) / c(3))) + c(4);

    % 2. Define a dense numerical search space centered around mu (coeffs(2))
    % A span of +/- 10*sigma guarantees capturing both half-max points even with high skew
    x_fine = linspace(coeffs(2) - 10*abs(sigma), coeffs(2) + 10*abs(sigma), 100000);
    y_fine = formSkewed_local(coeffs, x_fine);

    % 3. Isolate the pure baseline-subtracted peak profile
    baseline = coeffs(4);
    y_peak = y_fine - baseline; 

    % 4. Find the peak maximum and its matching half-maximum level
    % Works dynamically for both upright and inverted peaks via absolute value
    [peak_val, idx_peak] = max(abs(y_peak));
    half_max = peak_val / 2;

    % 5. Numerically cross-reference indices for the left and right split halves
    [~, idx_left] = min(abs(abs(y_peak(1:idx_peak)) - half_max));
    [~, idx_right] = min(abs(abs(y_peak(idx_peak:end)) - half_max));
    idx_right = idx_right + idx_peak - 1; % Account for index offset shift

    % 6. Map back to x values to compute the absolute geometric width
    fwhm_empirical = x_fine(idx_right) - x_fine(idx_left);

    % ---------------------------------------------------------
    % Center of Mass (True Mean) and Error Propagation
    % ---------------------------------------------------------
    mu = coeffs(2);
    mu_err = err(2);
    alpha = coeffs(5);
    alpha_err = err(5);

    delta = alpha / sqrt(1 + alpha^2);
    com = mu + sigma * delta * sqrt(2/pi);

    d_dmu = 1;
    d_dsigma = delta * sqrt(2/pi);
    d_dalpha = sigma * sqrt(2/pi) * (1 / (1 + alpha^2)^(3/2));
    com_err = sqrt( (d_dmu * mu_err)^2 + (d_dsigma * sigma_err)^2 + (d_dalpha * alpha_err)^2 );

    % ---------------------------------------------------------
    % Integral of the skewed Gaussian
    % ---------------------------------------------------------
    lineintegral = coeffs(1) * sigma * sqrt(2*pi);
    lineintegral_err = abs(lineintegral) * sqrt( ...
        (err(1)/coeffs(1))^2 + ...
        (sigma_err/sigma)^2 );

    % ---------------------------------------------------------
    % Update parameters for annotation
    % ---------------------------------------------------------

    % Replace mu with Center of Mass
    coeffs(2) = com;
    err(2) = com_err;

    % Replace sigma with FWHM for annotation
    coeffs(3) = fwhm_sym;
    err(3) = fwhm_sym_err;

    % ---------------------------------------------------------
    % Names and units
    % ---------------------------------------------------------

    if isempty(coeffNames)
        for i = 1:numel(coeffs)
            coeffNames{i} = ['Coeff ', num2str(i)];
        end
    end

    coeffNames{2} = 'Center of Mass';

    if isempty(coeffUnits)
        for i = 1:numel(coeffs)
            coeffUnits{i} = '';
        end
    end

    % ---------------------------------------------------------
    % Build annotation strings
    % ---------------------------------------------------------

    strs = cell(numel(coeffs)+2,1); % Expanded size for extra empirical FWHM line

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

    % Add line integral
    strs{end-1} = [ ...
        newline, ...
        'LineIntegral: ', ...
        unc_string(lineintegral,lineintegral_err), newline];

    % Add empirical FWHM line using the spatial x-axis units
    x_units = coeffUnits{2};
    strs{end} = ['Empirical FWHM: ', num2str(fwhm_empirical, '%.4g'), ' ', x_units];

    % ---------------------------------------------------------
    % Annotation
    % ---------------------------------------------------------

    an = annotation('textbox', dim, ...
        'String', strjoin(strs), ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white');

end

function an = myAnnotateDoubleSkewed(coeffs,err,coeffNames,coeffUnits)

    dim = [.6 .5 .3 .3];

    % =========================================================
    % coeffs =
    %
    % [A1,A2,sigma1,sigma2,Delta,x0,alpha1,alpha2,offset]
    %
    % Delta = separation between the Gaussian CENTERS
    % Peak-to-peak distance is calculated separately below.
    % =========================================================


    % ---------------------------------------------------------
    % Save original sigma values for peak calculation
    % ---------------------------------------------------------

    sigma1 = coeffs(3);
    sigma2 = coeffs(4);

    % ---------------------------------------------------------
    % Convert sigma -> FWHM for display
    % ---------------------------------------------------------

    coeffs(3) = coeffs(3)*2*sqrt(2*log(2));
    coeffs(4) = coeffs(4)*2*sqrt(2*log(2));

    err(3) = err(3)*2*sqrt(2*log(2));
    err(4) = err(4)*2*sqrt(2*log(2));


    % ---------------------------------------------------------
    % Calculate Gaussian areas
    %
    % Area = A*sigma*sqrt(2*pi)
    % ---------------------------------------------------------

    area1 = coeffs(1)*sigma1*sqrt(2*pi);
    area2 = coeffs(2)*sigma2*sqrt(2*pi);


    % ---------------------------------------------------------
    % Area ratio
    % ---------------------------------------------------------

    totalArea = abs(area1)+abs(area2);

    if totalArea ~= 0
        areaRatio = abs(area1)/totalArea;
    else
        areaRatio = NaN;
    end


    % ---------------------------------------------------------
    % Area uncertainties
    % ---------------------------------------------------------

    if coeffs(1) ~= 0

        area1_err = abs(area1)*sqrt( ...
            (err(1)/coeffs(1))^2 + ...
            (err(3)/coeffs(3))^2);

    else

        area1_err = NaN;

    end


    if coeffs(2) ~= 0

        area2_err = abs(area2)*sqrt( ...
            (err(2)/coeffs(2))^2 + ...
            (err(4)/coeffs(4))^2);

    else

        area2_err = NaN;

    end


    % ---------------------------------------------------------
    % Approximate uncertainty on area ratio
    % ---------------------------------------------------------

    if totalArea ~= 0 && ...
            ~isnan(area1_err) && ~isnan(area2_err)

        areaRatio_err = areaRatio * sqrt( ...
            (area1_err/abs(area1))^2 + ...
            ((area1_err+area2_err)/totalArea) );

    else

        areaRatio_err = NaN;

    end


    % =========================================================
    % ACTUAL PEAK POSITIONS
    %
    % For a skewed Gaussian:
    %
    % f(x) = exp(-(x-mu)^2/(2*sigma^2)) ...
    %        .* 2*normcdf(alpha*(x-mu)/sigma)
    %
    % The maximum is not generally at x = mu.
    %
    % We numerically find the maximum of each component.
    % =========================================================

    G1 = @(x) coeffs(1) .* ...
        exp(-(x-coeffs(6)).^2 ./ (2*coeffs(3)^2)) .* ...
        2 .* normcdf(coeffs(7).*(x-coeffs(6))./coeffs(3));
    
    G2 = @(x) coeffs(2) .* ...
        exp(-(x-(coeffs(6)+coeffs(5))).^2 ./ (2*coeffs(4)^2)) .* ...
        2 .* normcdf(coeffs(8).*(x-(coeffs(6)+coeffs(5)))./coeffs(4));
    
    xPeak1 = fminsearch(@(x) -G1(x), coeffs(6));
    xPeak2 = fminsearch(@(x) -G2(x), coeffs(6)+coeffs(5));

    yPeak1 = G1(xPeak1);
    yPeak2 = G2(xPeak2);




    % ---------------------------------------------------------
    % Peak-to-peak distance
    % ---------------------------------------------------------

    peakDistance = abs(xPeak2-xPeak1);


    % =========================================================
    % UNCERTAINTY ON PEAK-TO-PEAK DISTANCE
    %
    % Numerical propagation using the supplied coefficient
    % uncertainties.
    %
    % Only parameters that affect peak locations are relevant:
    %
    % sigma1  -> coeff 3
    % sigma2  -> coeff 4
    % Delta   -> coeff 5
    % x0      -> coeff 6
    % alpha1  -> coeff 7
    % alpha2  -> coeff 8
    % =========================================================

    % Function returning peak-to-peak distance from the
    % ORIGINAL parameterization.

    peakDistanceFunction = @(c) localPeakDistance(c);


    relevant = [3 4 5 6 7 8];

    variance = 0;

    for k = 1:numel(relevant)

        i = relevant(k);

        if err(i) == 0
            continue
        end

        % Work with original sigma uncertainties.
        % err(3) and err(4) were converted to FWHM above,
        % so convert them back here.

        if i == 3 || i == 4
            parameterErr = err(i)/(2*sqrt(2*log(2)));
        else
            parameterErr = err(i);
        end

        % Numerical derivative step

        step = max(abs(parameterErr)/10, ...
                   abs(coeffs(i))*1e-6);

        if step == 0
            step = 1e-6;
        end

        cPlus = coeffs;
        cMinus = coeffs;

        % coeffs(3:4) currently contain FWHM, so convert them
        % back to sigma before passing to localPeakDistance.

        cPlus(3) = sigma1;
        cPlus(4) = sigma2;

        cMinus(3) = sigma1;
        cMinus(4) = sigma2;

        cPlus(i) = cPlus(i) + step;
        cMinus(i) = cMinus(i) - step;

        dPlus = localPeakDistance(cPlus);
        dMinus = localPeakDistance(cMinus);

        derivative = (dPlus-dMinus)/(2*step);

        variance = variance + ...
            (derivative*parameterErr)^2;

    end

    peakDistance_err = sqrt(variance);


    % ---------------------------------------------------------
    % Default names and units
    % ---------------------------------------------------------

    if isempty(coeffNames)

        for i = 1:numel(coeffs)

            coeffNames{i} = ...
                ['Coeff ',num2str(i)];

        end

    end


    if isempty(coeffUnits)

        for i = 1:numel(coeffs)

            coeffUnits{i} = '';

        end

    end


    % =========================================================
    % Build annotation
    % =========================================================

    strs = cell(numel(coeffs)+1,1);


    for i = 1:numel(coeffs)

        if i < numel(coeffs)

            strs{i} = [ ...
                coeffNames{i},': ',...
                unc_string(coeffs(i),err(i)),...
                ' ',coeffUnits{i},newline];

        else

            strs{i} = [ ...
                coeffNames{i},': ',...
                unc_string(coeffs(i),err(i)),...
                ' ',coeffUnits{i}];

        end

    end


    % ---------------------------------------------------------
    % Add peak positions
    % ---------------------------------------------------------

    % strs{numel(coeffs)+1} = [ ...
    %     newline,...
    %     'Peak 1: ',...
    %     num2str(yPeak1)];
% 
    % strs{numel(coeffs)+2} = [ ...
    %     'Peak 2: ',...
    %     num2str(yPeak2)];


    % ---------------------------------------------------------
    % Add peak-to-peak distance
    % ---------------------------------------------------------

    strs{numel(coeffs)+1} = [ ...
        'Peak-to-Peak Distance: ',...
        unc_string(peakDistance,peakDistance_err)];


    % ---------------------------------------------------------
    % Add areas
    % ---------------------------------------------------------

    % strs{numel(coeffs)+4} = [ ...
    %     newline,...
    %     'Area 1: ',...
    %     num2str(area1),newline,...
    %     'Area 2: ',...
    %     num2str(area2),newline,...
    %     'Area Ratio: ',...
    %     unc_string(areaRatio,areaRatio_err)];


    % ---------------------------------------------------------
    % Annotation
    % ---------------------------------------------------------

    an = annotation('textbox',dim,...
        'String',strjoin(strs),...
        'FitBoxToText','on',...
        'BackgroundColor','white');


end


% =============================================================
% LOCAL FUNCTION
% Calculate actual peak-to-peak distance.
%
% Input c uses the ORIGINAL parameterization:
%
% [A1,A2,sigma1,sigma2,Delta,x0,alpha1,alpha2,offset]
% =============================================================

function d = localPeakDistance(c)

    sigma1 = c(3);
    sigma2 = c(4);

    mu1 = c(6);
    mu2 = c(6) + c(5);

    alpha1 = c(7);
    alpha2 = c(8);


    % Skewed Gaussian #1

    f1 = @(x) ...
        exp(-(x-mu1).^2 ./ (2*sigma1^2)) .* ...
        2 .* normcdf(alpha1 .* (x-mu1) ./ sigma1);


    % Skewed Gaussian #2

    f2 = @(x) ...
        exp(-(x-mu2).^2 ./ (2*sigma2^2)) .* ...
        2 .* normcdf(alpha2 .* (x-mu2) ./ sigma2);


    % Search ranges

    range1 = 10*sigma1;
    range2 = 10*sigma2;


    % Find actual maxima

    [~,peak1] = fminbnd( ...
        @(x)-f1(x), ...
        mu1-range1, ...
        mu1+range1);


    [~,peak2] = fminbnd( ...
        @(x)-f2(x), ...
        mu2-range2, ...
        mu2+range2);


    % Peak-to-peak distance

    d = abs(peak2-peak1);

end

