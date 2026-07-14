function indivDataset = param_extract_sfi_doubleSkew_TwoLib(indivDataset, leftTransitionName, rightTransitionName, libFile)
% Double-peak fit using TWO LIBRARY templates on counts/bin in the voltage domain.
%
% - LEFT peak:  [muL, sigmaL, fL, tauL] fixed from library; fit only Aleft
% - RIGHT peak: [muR, sigmaR, fR, tauR] fixed from library; fit only Aright
% - Shared baseline b0 (counts/bin)
%
% Parameter vector q:
%   q = [Aleft, Aright, b0]
%
% Output fields (per scan i):
%   indiv.sfiIntegralLeft{j}      = Aleft   (counts above baseline)
%   indiv.sfiIntegralRight{j}     = Aright  (counts above baseline)
%   indiv.sfiFitA_left(j), indiv.sfiFitA_right(j), indiv.sfiFitB0(j)
%   indiv.sfiDoubleFitInitParams{j} = [Aleft0 Aright0 b0_0]
%   indiv.sfiDoubleFitParams{j}     = [Aleft Aright b0]
%   indiv.sfiDoubleFitCounts{j}     = yhat total
%   indiv.sfiDoubleFitCountsLeft{j} = left component (no baseline)
%   indiv.sfiDoubleFitCountsRight{j}= right component (no baseline)

    % ------------------------------------------------------------
    % User option: enable summing SFI by voltage for low SFI counts
    % ------------------------------------------------------------
    sumSFIByVoltageForLowSFICounts_bool = false;

    d = dialog( ...
        'Name', 'SFI Fit Options', ...
        'Position', [500 500 360 140], ...
        'WindowStyle', 'modal');

    cb = uicontrol( ...
        'Parent', d, ...
        'Style', 'checkbox', ...
        'String', 'Enable Sum SFI by Voltage for low SFI counts', ...
        'Value', 0, ...
        'Position', [20 75 320 30], ...
        'HorizontalAlignment', 'left');

    uicontrol( ...
        'Parent', d, ...
        'Style', 'pushbutton', ...
        'String', 'OK', ...
        'Position', [135 20 90 30], ...
        'Callback', @(src,evt) uiresume(d));

    uiwait(d);

    if isvalid(d)
        sumSFIByVoltageForLowSFICounts_bool = logical(cb.Value);
        delete(d);
    end
    if nargin < 4 || isempty(libFile)
        writerPath = which('saveSFIParamLibrary_UI');
        if ~isempty(writerPath)
            libFile = fullfile(fileparts(writerPath), 'sfi_fit_parameter_library.txt');
        else
            here = fileparts(mfilename('fullpath'));
            libFile = fullfile(here, 'sfi_fit_parameter_library.txt');
        end
    end
    % ---- Standardized plot-data export controls ----
    savePlotDataBool = true;
    plotDataExportFolder = "E:\SFI_Exports";
    dataSetType = "DoubleFitLibraryPlusLibrary_SummedByVoltage";

    % ---- Load both templates once ----
    [muL, sigmaL, fL, tauL] = local_get_template_doubleSkew(libFile, leftTransitionName);
    [muR, sigmaR, fR, tauR] = local_get_template_doubleSkew(libFile, rightTransitionName);

    % Safety / convention check
    if ~(isfinite(muL) && isfinite(muR))
        error('Non-finite library template center(s).');
    end

    if muL >= muR
        warning(['Left template center is not less than right template center. ', ...
                 'Proceeding anyway, but check that the selected transitions are assigned correctly.']);
    end
    % =========================================================
    % Optional branch: sum SFI spectra by repeated imagevcoAtom values
    % =========================================================
    if sumSFIByVoltageForLowSFICounts_bool

        T = [];      % unique imagevcoAtom values
        X = {};      % averaged SFI spectra corresponding to each T entry
        Xsum = {};   % raw summed SFI spectra before averaging
        Nsum = [];   % number of spectra contributing to each T entry
        
        debugPlotSummedSFI_bool = true;
        for i = 1:numel(indivDataset)

            indiv = indivDataset{i};
            if isempty(indiv) || ~isstruct(indiv), continue; end
            if ~isfield(indiv, 'CounterMCS') || indiv.CounterMCS < 1, continue; end
            if ~isfield(indiv, 'mcsSpectra') || isempty(indiv.mcsSpectra), continue; end
            if ~isfield(indiv, 'imagevcoAtom') || isempty(indiv.imagevcoAtom), continue; end

            Nshots = min([indiv.CounterMCS, numel(indiv.mcsSpectra), numel(indiv.imagevcoAtom)]);

            for j = 1:Nshots
                if isempty(indiv.mcsSpectra{j})
                    continue;
                end

                spec = indiv.mcsSpectra{j};
                if size(spec,2) < 2
                    continue;
                end

                thisT = indiv.imagevcoAtom(j);

                if ~isfinite(thisT)
                    continue;
                end

                % keep the full SFI array for now
                thisX = spec;

                % search for matching time already in T
                n = find(T == thisT, 1, 'first');

                if isempty(n)
                    T(end+1,1) = thisT;       %#ok<AGROW>
                
                    Xsum{end+1,1} = thisX;    %#ok<AGROW>
                    Nsum(end+1,1) = 1;        %#ok<AGROW>
                
                else
                    if isequal(size(Xsum{n}), size(thisX)) && all(Xsum{n}(:,1) == thisX(:,1))
                        Xsum{n}(:,2) = Xsum{n}(:,2) + thisX(:,2);
                        Nsum(n) = Nsum(n) + 1;
                    else
                        warning(['Skipping average for scan %d shot %d: SFI spectrum/time axis ', ...
                                 'does not match existing grouped spectrum for imagevcoAtom = %.9g'], ...
                                 i, j, thisT);
                    end
                end
            end
        end

        % sort by T for cleanliness
        [T, ordT] = sort(T);
        Xsum = Xsum(ordT);
        Nsum = Nsum(ordT);
        % Convert summed spectra into averaged spectra
        X = cell(size(Xsum));
        
        for n = 1:numel(Xsum)
            X{n} = Xsum{n};
            if Nsum(n) > 0
                X{n}(:,2) = Xsum{n}(:,2) ./ Nsum(n);
            else
                X{n}(:,2) = NaN(size(Xsum{n}(:,2)));
            end
        end
        % ---------------------------------------------------------
        % Debug plots for grouped / summed SFI spectra
        % ---------------------------------------------------------
        if debugPlotSummedSFI_bool

            nSpec = numel(T);
            nCols = ceil(sqrt(nSpec));
            nRows = ceil(nSpec / nCols);
        
            figure('Name', 'Averaged grouped SFI debug (all T)', 'Color', 'w');
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        
            for n = 1:nSpec
        
                specPlot = X{n};
                if isempty(specPlot) || size(specPlot,2) < 2
                    continue;
                end
        
                nexttile;
                hold on; box on;
        
                plot(specPlot(:,1), specPlot(:,2), 'k.-', ...
                     'LineWidth', 1.0, 'MarkerSize', 6);
        
                title(sprintf('T = %.6g', T(n)));
                xlabel('Time');
                ylabel('Counts');
        
                set(gca, 'FontSize', 8);
        
                hold off;
            end
        end
        % =========================================================
        % Fit all grouped / summed SFI spectra in X vs T
        % using the same strategy as the non-summed branch
        % =========================================================
        Ngroup = numel(T);

        groupedFitInitParams   = cell(Ngroup,1);
        groupedFitParams       = cell(Ngroup,1);
        groupedFitCounts       = cell(Ngroup,1);
        groupedFitCountsLeft   = cell(Ngroup,1);
        groupedFitCountsRight  = cell(Ngroup,1);

        groupedIntegralLeft    = nan(Ngroup,1);
        groupedIntegralRight   = nan(Ngroup,1);
        groupedIntegralTotal   = nan(Ngroup,1);

        groupedIntegralLeft_err  = nan(Ngroup,1);
        groupedIntegralRight_err = nan(Ngroup,1);
        groupedB0_err            = nan(Ngroup,1);

        groupedChi2            = nan(Ngroup,1);
        groupedChi2Red         = nan(Ngroup,1);
        groupedFitValid        = false(Ngroup,1);

        for n = 1:Ngroup

            spec = X{n};
            if isempty(spec) || size(spec,2) < 2
                continue;
            end

            t = spec(:,1);
            y = spec(:,2);   % averaged counts/bin
            
            nContrib = Nsum(n);

            if numel(t) < 5
                continue;
            end

            t = t(:);
            y = y(:);

            % ---- convert to V axis ----
            [V, dV] = timeToVoltageAxis(t);

            % enforce increasing ordering
            [V, ord] = sort(V);
            y  = y(ord);
            dV = dV(ord); %#ok<NASGU>
            
            % Uncertainty of averaged counts/bin:
            % Y_i = summed counts/bin before averaging
            % N   = number of spectra averaged
            % sigma_yi = sqrt(max(Y_i, 1)) / N
            if isfinite(nContrib) && nContrib > 0
                Ysum_i = y .* nContrib;
                ySigma = sqrt(max(Ysum_i, 1)) ./ nContrib;
            else
                ySigma = ones(size(y));
            end
            
            if numel(V) < 5, continue; end
            if any(~isfinite(V)) || any(~isfinite(y)), continue; end

            % ---- baseline guess ----
            b0_0 = median(y(1:min(20, numel(y))));

            % counts above baseline for initial guesses
            yb = y - b0_0;
            yb(yb < 0) = 0;

            % ---- unit-amplitude library models in counts/bin ----
            mL = skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], V);
            mR = skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], V);

            if any(~isfinite(mL)) || any(~isfinite(mR))
                continue;
            end

            % ---- initial amplitude guess from linear least squares ----
            X0 = [mL(:), mR(:)];
            c0 = X0 \ yb;
            if numel(c0) < 2 || any(~isfinite(c0))
                c0 = [0; 0];
            end

            Aleft_0  = max(c0(1), 0);
            Aright_0 = max(c0(2), 0);

            q0 = [Aleft_0, Aright_0, b0_0];
            groupedFitInitParams{n} = q0(:).';

            % rough amplitude cap scale
            Atot = sum(yb);
            if ~isfinite(Atot) || Atot <= 0
                Atot = sum(max(y - b0_0, 0));
            end
            Atot = max(Atot, 1);
            Amax = 5 * Atot;

            lb = [0,    0,    -Inf];
            ub = [Amax, Amax,  Inf];

            model = @(q, Vaxis) ( ...
                  q(1) .* skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], Vaxis) ...
                + q(2) .* skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], Vaxis) ...
                + q(3) );

            opt = optimoptions('lsqcurvefit', 'Display', 'off');

            try
                [qfit, resnorm, residual, exitflag, output, lambda, J] = ...
                    lsqcurvefit(model, q0, V, y, lb, ub, opt); %#ok<ASGLU>
            catch
                continue;
            end

            % ---- evaluate and store ----
            yhat_left  = qfit(1) .* skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], V);
            yhat_right = qfit(2) .* skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], V);
            yhat_tot   = yhat_left + yhat_right + qfit(3);

            p = 3;
            N = numel(y);

            resid = y - yhat_tot;
            
            validChiPoints = isfinite(resid) & isfinite(ySigma) & (ySigma > 0);
            
            chi2 = sum((resid(validChiPoints) ./ ySigma(validChiPoints)).^2);
            
            dof = nnz(validChiPoints) - p;
            if dof > 0
                chi2red = chi2 / dof;
            else
                chi2red = NaN;
            end
                        
            fprintf('SFI grouped fit diagnostic: T = %.6g | Nsum = %d | RMS resid = %.4g | mean ySigma = %.4g | chi2red = %.4g\n', ...
                T(n), nContrib, sqrt(mean(resid(validChiPoints).^2)), mean(ySigma(validChiPoints)), chi2red);

            groupedChi2(n)        = chi2;
            groupedChi2Red(n)     = chi2red;
            groupedIntegralLeft(n)  = qfit(1);
            groupedIntegralRight(n) = qfit(2);
            groupedIntegralTotal(n) = qfit(1) + qfit(2);

            % ---- grouped parameter uncertainties from covariance matrix ----
            qerr = [NaN; NaN; NaN];
            
            if ~isempty(J)
                dof_cov = N - numel(qfit);
                if dof_cov > 0
                    JTJ = full(J' * J);
                    if all(isfinite(JTJ(:)))
                        sigma2 = resnorm / dof_cov;
                    
                        if rcond(JTJ) > 1e-12
                            cov_q = sigma2 * inv(JTJ);
                        else
                            cov_q = sigma2 * pinv(JTJ);
                        end
                    
                        if all(isfinite(cov_q(:)))
                            qerr = sqrt(diag(cov_q));
                        end
                    end
                end
            end
            
            groupedIntegralLeft_err(n)  = qerr(1);
            groupedIntegralRight_err(n) = qerr(2);
            groupedB0_err(n)            = qerr(3);

            groupedFitParams{n}      = qfit(:).';
            groupedFitCounts{n}      = yhat_tot;
            groupedFitCountsLeft{n}  = yhat_left;
            groupedFitCountsRight{n} = yhat_right;
            groupedFitValid(n)       = true;
        end
        % =========================================================
        % Grouped ROI-style outputs based on grouped Aleft / Aright
        % =========================================================
        grouped_roi1 = groupedIntegralLeft;
        grouped_roi2 = groupedIntegralRight;

        grouped_roi1_ratio = nan(size(T));
        grouped_roi2_ratio = nan(size(T));

        grouped_roi1_ratio_err = nan(size(T));
        grouped_roi2_ratio_err = nan(size(T));

        denomLR = grouped_roi1 + grouped_roi2;

        validROI = groupedFitValid & ...
                   isfinite(grouped_roi1) & ...
                   isfinite(grouped_roi2) & ...
                   isfinite(groupedIntegralLeft_err) & ...
                   isfinite(groupedIntegralRight_err) & ...
                   isfinite(denomLR) & ...
                   (denomLR > 0);

        grouped_roi1_ratio(validROI) = grouped_roi1(validROI) ./ denomLR(validROI);
        grouped_roi2_ratio(validROI) = grouped_roi2(validROI) ./ denomLR(validROI);

        % Error propagation for:
        % f1 = Aleft / (Aleft + Aright)
        % f2 = Aright / (Aleft + Aright)
        %
        % ignoring covariance between Aleft and Aright for now

        A = grouped_roi1(validROI);
        B = grouped_roi2(validROI);
        sA = groupedIntegralLeft_err(validROI);
        sB = groupedIntegralRight_err(validROI);
        D = A + B;

        df1_dA =  B ./ (D.^2);
        df1_dB = -A ./ (D.^2);

        df2_dA = -B ./ (D.^2);
        df2_dB =  A ./ (D.^2);

        grouped_roi1_ratio_err(validROI) = sqrt( (df1_dA.^2).*(sA.^2) + ...
                                                 (df1_dB.^2).*(sB.^2) );

        grouped_roi2_ratio_err(validROI) = sqrt( (df2_dA.^2).*(sA.^2) + ...
                                                 (df2_dB.^2).*(sB.^2) );
                if savePlotDataBool
            export_two_library_summed_by_voltage_csv( ...
                plotDataExportFolder, ...
                dataSetType, ...
                leftTransitionName, ...
                rightTransitionName, ...
                libFile, ...
                T, X, Xsum, Nsum, ...
                groupedFitInitParams, ...
                groupedFitParams, ...
                groupedFitCounts, ...
                groupedFitCountsLeft, ...
                groupedFitCountsRight, ...
                groupedIntegralLeft, ...
                groupedIntegralRight, ...
                groupedIntegralTotal, ...
                groupedIntegralLeft_err, ...
                groupedIntegralRight_err, ...
                groupedB0_err, ...
                groupedChi2, ...
                groupedChi2Red, ...
                groupedFitValid, ...
                grouped_roi1_ratio, ...
                grouped_roi2_ratio, ...
                grouped_roi1_ratio_err, ...
                grouped_roi2_ratio_err, ...
                [muL, sigmaL, fL, tauL], ...
                [muR, sigmaR, fR, tauR]);
        end
        % =========================================================
        % Plot 1: all grouped fits in one tiled figure
        % =========================================================
        plotIdx = find(groupedFitValid);
        if ~isempty(plotIdx)

            nValid = numel(plotIdx);
            nCols = ceil(sqrt(nValid));
            nRows = ceil(nValid / nCols);

            figure('Name', sprintf('Grouped SFI fits (%d grouped spectra)', nValid), ...
                   'Color', 'w');
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

            for kk = 1:nValid
                n = plotIdx(kk);

            spec = X{n};
            t_plot = spec(:,1);
            y_plot = spec(:,2);
            
            nContrib_plot = Nsum(n);
            
            [V_plot, dV_plot] = timeToVoltageAxis(t_plot); %#ok<NASGU>
            [V_plot, ord] = sort(V_plot);
            y_plot = y_plot(ord);
            
            % Uncertainty of averaged counts/bin
            if isfinite(nContrib_plot) && nContrib_plot > 0
                yerr_plot = sqrt(max(y_plot .* nContrib_plot, 1)) ./ nContrib_plot;
            else
                yerr_plot = ones(size(y_plot));
            end

                yhat_left  = groupedFitCountsLeft{n};
                yhat_right = groupedFitCountsRight{n};
                yhat_tot   = groupedFitCounts{n};

                nexttile;
                hold on; box on;

                plot(V_plot, yhat_tot,   'r-',  'LineWidth', 1.5, 'DisplayName', 'Total fit');
                plot(V_plot, yhat_left,  '--',  'LineWidth', 1.0, 'DisplayName', 'Left lib');
                plot(V_plot, yhat_right, '--',  'LineWidth', 1.0, 'DisplayName', 'Right lib');
                errorbar(V_plot, y_plot, yerr_plot, ...
                    'ko', ...
                    'MarkerSize', 2, ...
                    'MarkerFaceColor', 'k', ...
                    'LineWidth', 0.8, ...
                    'CapSize', 3, ...
                    'DisplayName', 'Data');

                title(sprintf('T = %.6g | \\chi^2_\\nu = %.3f', T(n), groupedChi2Red(n)));
                xlabel('Voltage (V)');
                ylabel('Counts/bin');
                set(gca, 'FontSize', 8);

                if kk == 1
                    legend('Location', 'best');
                end

                hold off;
            end

            validChi = isfinite(groupedChi2Red(plotIdx));
            if any(validChi)
                meanChi = mean(groupedChi2Red(plotIdx(validChi)));
                sgtitle(sprintf('Grouped summed SFI fits | mean \\chi^2_\\nu = %.3f', meanChi), ...
                        'FontSize', 14);
            else
                sgtitle('Grouped summed SFI fits', 'FontSize', 14);
            end
        end
                % =========================================================
        % Debug single grouped SFI fit plot
        % =========================================================
        if debugPlotSummedSFI_bool

            debugSingleIdx = 12;   % choose which grouped spectrum to inspect

            if debugSingleIdx <= numel(T) && groupedFitValid(debugSingleIdx)

                n = debugSingleIdx;

                spec = X{n};
                t_plot = spec(:,1);
                y_plot = spec(:,2);

                nContrib_plot = Nsum(n);

                [V_plot, dV_plot] = timeToVoltageAxis(t_plot); %#ok<NASGU>
                [V_plot, ord] = sort(V_plot);
                y_plot = y_plot(ord);

                % Uncertainty of averaged counts/bin:
                % Y_i = summed counts/bin before averaging
                % sigma_yi = sqrt(max(Y_i, 1)) / N
                if isfinite(nContrib_plot) && nContrib_plot > 0
                    Ysum_plot = y_plot .* nContrib_plot;
                    yerr_plot = sqrt(max(Ysum_plot, 1)) ./ nContrib_plot;
                else
                    yerr_plot = ones(size(y_plot));
                end

                yhat_left  = groupedFitCountsLeft{n};
                yhat_right = groupedFitCountsRight{n};
                yhat_tot   = groupedFitCounts{n};

                Aleft  = groupedIntegralLeft(n);
                Aright = groupedIntegralRight(n);
                Atotal = groupedIntegralTotal(n);

                AleftErr  = groupedIntegralLeft_err(n);
                ArightErr = groupedIntegralRight_err(n);

                chi2    = groupedChi2(n);
                chi2red = groupedChi2Red(n);

                ratio1 = grouped_roi1_ratio(n);
                ratio2 = grouped_roi2_ratio(n);

                ratio1Err = grouped_roi1_ratio_err(n);
                ratio2Err = grouped_roi2_ratio_err(n);

                figSingle = figure('Name', sprintf('Debug grouped SFI fit n=%d', n), ...
                                   'Color', 'w');
                hold on; box on;

                errorbar(V_plot, y_plot, yerr_plot, ...
                    'ko', ...
                    'MarkerSize', 4, ...
                    'MarkerFaceColor', 'k', ...
                    'LineWidth', 1.0, ...
                    'CapSize', 5, ...
                    'DisplayName', 'Data');

                plot(V_plot, yhat_tot, ...
                    'r-', ...
                    'LineWidth', 2.0, ...
                    'DisplayName', 'Total fit');

                plot(V_plot, yhat_left, ...
                    '--', ...
                    'LineWidth', 1.5, ...
                    'DisplayName', 'Left lib');

                plot(V_plot, yhat_right, ...
                    '--', ...
                    'LineWidth', 1.5, ...
                    'DisplayName', 'Right lib');

                xlabel('Voltage (V)');
                ylabel('Averaged counts/bin');
                title(sprintf('Grouped SFI fit debug | n = %d | T = %.6g | Nsum = %d', ...
                    n, T(n), nContrib_plot));

                legend('Location', 'best');
                set(gca, 'FontSize', 14);

                statText = sprintf([ ...
                    'n = %d\n', ...
                    'T = %.6g\n', ...
                    'Nsum = %d\n', ...
                    '\\chi^2 = %.4g\n', ...
                    '\\chi^2_\\nu = %.4g\n', ...
                    'A_L = %.4g +/- %.4g\n', ...
                    'A_R = %.4g +/- %.4g\n', ...
                    'A_L + A_R = %.4g\n', ...
                    'ROI1 = %.4g +/- %.4g\n', ...
                    'ROI2 = %.4g +/- %.4g'], ...
                    n, T(n), nContrib_plot, ...
                    chi2, chi2red, ...
                    Aleft, AleftErr, ...
                    Aright, ArightErr, ...
                    Atotal, ...
                    ratio1, ratio1Err, ...
                    ratio2, ratio2Err);

                annotation(figSingle, 'textbox', [0.62 0.55 0.30 0.35], ...
                    'String', statText, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'white', ...
                    'EdgeColor', 'black', ...
                    'FontSize', 11);

                hold off;

            else
                warning('Debug single grouped SFI plot skipped: index %d is invalid or not fit-valid.', debugSingleIdx);
            end
        end

        % =========================================================
        % Plot 2: peak1 / (peak1 + peak2) vs T
        % =========================================================
        ratio12 = nan(size(T));
        denom12 = groupedIntegralLeft + groupedIntegralRight;

        validRatio = groupedFitValid & isfinite(groupedIntegralLeft) & ...
                     isfinite(groupedIntegralRight) & isfinite(denom12) & ...
                     (denom12 > 0);

        ratio12(validRatio) = groupedIntegralLeft(validRatio) ./ denom12(validRatio);

        if any(validRatio)
            figure('Name', 'Grouped peak fraction vs T', 'Color', 'w');
            hold on; box on;

            plot(T(validRatio), ratio12(validRatio), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);

            xlabel('T');
            ylabel('Peak1 / (Peak1 + Peak2)');
            title('Grouped averaged SFI: Peak1 fraction vs T');
            ylim([0 1]);
            set(gca, 'FontSize', 14);

            hold off;
        end

                % store grouped result and grouped fit outputs
        for i = 1:numel(indivDataset)
            if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
                continue;
            end

            indivDataset{i}.sumSFIByVoltageForLowSFICounts_bool = true;

            indivDataset{i}.summedSFIByVoltage_T = T;
            indivDataset{i}.summedSFIByVoltage_X = X;
            indivDataset{i}.summedSFIByVoltage_Xsum = Xsum;
            indivDataset{i}.summedSFIByVoltage_Nsum = Nsum;

            indivDataset{i}.summedSFIByVoltage_fitInitParams   = groupedFitInitParams;
            indivDataset{i}.summedSFIByVoltage_fitParams       = groupedFitParams;
            indivDataset{i}.summedSFIByVoltage_fitCounts       = groupedFitCounts;
            indivDataset{i}.summedSFIByVoltage_fitCountsLeft   = groupedFitCountsLeft;
            indivDataset{i}.summedSFIByVoltage_fitCountsRight  = groupedFitCountsRight;

            indivDataset{i}.summedSFIByVoltage_integralLeft    = groupedIntegralLeft;
            indivDataset{i}.summedSFIByVoltage_integralRight   = groupedIntegralRight;
            indivDataset{i}.summedSFIByVoltage_integralTotal   = groupedIntegralTotal;

            indivDataset{i}.summedSFIByVoltage_chi2            = groupedChi2;
            indivDataset{i}.summedSFIByVoltage_chi2Red         = groupedChi2Red;
            indivDataset{i}.summedSFIByVoltage_fitValid        = groupedFitValid;
            indivDataset{i}.summedSFIByVoltage_ratio12         = ratio12;

            indivDataset{i}.summedSFIByVoltage_roi1       = grouped_roi1;
            indivDataset{i}.summedSFIByVoltage_roi2       = grouped_roi2;
            indivDataset{i}.summedSFIByVoltage_roi1_ratio = grouped_roi1_ratio;
            indivDataset{i}.summedSFIByVoltage_roi2_ratio = grouped_roi2_ratio;
            
            % Compatibility fields for base_fit / lorentzian_lineshape
            indivDataset{i}.imagevcoAtom = T;
            
            indivDataset{i}.sfiIntegral_roi1       = grouped_roi1;
            indivDataset{i}.sfiIntegral_roi2       = grouped_roi2;
            indivDataset{i}.sfiIntegral_roi1_ratio = grouped_roi1_ratio;
            indivDataset{i}.sfiIntegral_roi2_ratio = grouped_roi2_ratio;

            indivDataset{i}.summedSFIByVoltage_integralLeft_err  = groupedIntegralLeft_err;
            indivDataset{i}.summedSFIByVoltage_integralRight_err = groupedIntegralRight_err;
            indivDataset{i}.summedSFIByVoltage_b0_err            = groupedB0_err;

            indivDataset{i}.summedSFIByVoltage_roi1_ratio_err = grouped_roi1_ratio_err;
            indivDataset{i}.summedSFIByVoltage_roi2_ratio_err = grouped_roi2_ratio_err;
        end

        fprintf('Constructed grouped SFI spectra: %d unique imagevcoAtom values found.\n', numel(T));

        return;
    end

    for i = 1:numel(indivDataset)

        indiv = indivDataset{i};
        if isempty(indiv) || ~isstruct(indiv), continue; end
        if ~isfield(indiv,'CounterMCS') || indiv.CounterMCS < 1, continue; end
        if ~isfield(indiv,'mcsSpectra') || isempty(indiv.mcsSpectra), continue; end

        if ~isfield(indiv,'imagevcoAtom') || isempty(indiv.imagevcoAtom)
            warning('Scan %d skipped: missing imagevcoAtom.', i);
            continue;
        end
        
        Nshots = min([indiv.CounterMCS, numel(indiv.mcsSpectra), numel(indiv.imagevcoAtom)]);
        
        % Compatibility for getxy/base_fit
        indiv.CounterAtom = Nshots;
        indiv.imagevcoAtom = indiv.imagevcoAtom(1:Nshots);
        indiv.imagevcoAtom = indiv.imagevcoAtom(:).';

        % ---- metadata ----
        indiv.sfiDoubleFit_mode            = 'TwoLib';
        indiv.sfiDoubleFit_leftTransition  = char(leftTransitionName);
        indiv.sfiDoubleFit_rightTransition = char(rightTransitionName);
        indiv.sfiDoubleFit_leftTemplate    = [muL, sigmaL, fL, tauL];
        indiv.sfiDoubleFit_rightTemplate   = [muR, sigmaR, fR, tauR];

        % ---- outputs ----
        % Newer fields
        indiv.sfiIntegral           = nan(1, Nshots);
        indiv.sfiIntegralLeft       = nan(1, Nshots);
        indiv.sfiIntegralRight      = nan(1, Nshots);

        indiv.sfiFitA_left          = nan(1, Nshots);
        indiv.sfiFitA_right         = nan(1, Nshots);
        indiv.sfiFitB0              = nan(1, Nshots);

        indiv.sfiDoubleFitInitParams  = cell(1, Nshots);
        indiv.sfiDoubleFitParams      = cell(1, Nshots);
        indiv.sfiDoubleFitCounts      = cell(1, Nshots);
        indiv.sfiDoubleFitCountsLeft  = cell(1, Nshots);
        indiv.sfiDoubleFitCountsRight = cell(1, Nshots);

        indiv.sfiDoubleFitChi2Red   = nan(1, Nshots);
        indiv.sfiDoubleFitChi2      = nan(1, Nshots);

        % Legacy fields expected by Rabi analysis / base_fit
        indiv.sfiIntegral_roi1       = nan(1, Nshots);
        indiv.sfiIntegral_roi2       = nan(1, Nshots);
        indiv.sfiIntegral_roi1_ratio = nan(1, Nshots);
        indiv.sfiIntegral_roi2_ratio = nan(1, Nshots);

        indiv.sfiDoubleFitValid      = false(Nshots,1);

        for j = 1:Nshots
            if isempty(indiv.mcsSpectra{j}), continue; end
            spec = indiv.mcsSpectra{j};
            if size(spec,2) < 2, continue; end

            t = spec(:,1);
            y = spec(:,2);   % counts/bin
            if numel(t) < 5, continue; end

            t = t(:);
            y = y(:);

            % ---- convert to V axis ----
            [V, dV] = timeToVoltageAxis(t);

            % enforce increasing ordering
            [V, ord] = sort(V);
            y  = y(ord);
            dV = dV(ord); %#ok<NASGU>

            if numel(V) < 5, continue; end
            if any(~isfinite(V)) || any(~isfinite(y)), continue; end

            % ---- baseline guess ----
            b0_0 = median(y(1:min(20,numel(y))));

            % counts above baseline for initial guesses
            yb = y - b0_0;
            yb(yb < 0) = 0;

            % ---- unit-amplitude library models in counts/bin ----
            mL = skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], V);
            mR = skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], V);

            if any(~isfinite(mL)) || any(~isfinite(mR))
                continue;
            end

            % ---- initial amplitude guess from linear least squares on baseline-subtracted data ----
            X0 = [mL(:), mR(:)];
            c0 = X0 \ yb;
            if numel(c0) < 2 || any(~isfinite(c0))
                c0 = [0; 0];
            end

            Aleft_0  = max(c0(1), 0);
            Aright_0 = max(c0(2), 0);

            q0 = [Aleft_0, Aright_0, b0_0];
            indiv.sfiDoubleFitInitParams{j} = q0(:).';

            % rough amplitude cap scale
            Atot = sum(yb);
            if ~isfinite(Atot) || Atot <= 0
                Atot = sum(max(y - b0_0, 0));
            end
            Atot = max(Atot, 1);
            Amax = 5 * Atot;

            lb = [0,    0,    -Inf];
            ub = [Amax, Amax,  Inf];

            % model with both fixed templates
            model = @(q, Vaxis) ( ...
                  q(1) .* skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], Vaxis) ...
                + q(2) .* skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], Vaxis) ...
                + q(3) );

            opt = optimoptions('lsqcurvefit', 'Display','off');

            try
                qfit = lsqcurvefit(model, q0, V, y, lb, ub, opt);
            catch
                continue;
            end

            % ---- evaluate & store ----
            yhat_left  = qfit(1) .* skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], V);
            yhat_right = qfit(2) .* skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], V);
            yhat       = yhat_left + yhat_right + qfit(3);

            % ---- reduced chi-squared (unit uncertainty) ----
            p = 3;                 % [Aleft, Aright, b0]
            N = numel(y);

            resid = y - yhat;
            chi2 = sum(resid.^2);

            dof = N - p;
            if dof > 0
                chi2red = chi2 / dof;
            else
                chi2red = NaN;
            end


            indiv.sfiDoubleFitChi2(j)    = chi2;
            indiv.sfiDoubleFitChi2Red(j) = chi2red;

            Aleft  = qfit(1);
            Aright = qfit(2);
            b0     = qfit(3);

            indiv.sfiIntegralLeft(j)   = Aleft;
            indiv.sfiIntegralRight(j)  = Aright;
            indiv.sfiIntegral(j)       = Aleft + Aright;

            % Legacy ROI-style variables for downstream Rabi analysis
            indiv.sfiIntegral_roi1(j) = Aleft;
            indiv.sfiIntegral_roi2(j) = Aright;

            AtotLR = Aleft + Aright;
            if isfinite(AtotLR) && AtotLR > 0
                indiv.sfiIntegral_roi1_ratio(j) = Aleft  / AtotLR;
                indiv.sfiIntegral_roi2_ratio(j) = Aright / AtotLR;
            else
                indiv.sfiIntegral_roi1_ratio(j) = NaN;
                indiv.sfiIntegral_roi2_ratio(j) = NaN;
            end

            indiv.sfiFitA_left(j)      = Aleft;
            indiv.sfiFitA_right(j)     = Aright;
            indiv.sfiFitB0(j)          = b0;

            indiv.sfiDoubleFitParams{j}      = qfit(:).';
            indiv.sfiDoubleFitCounts{j}      = yhat;
            indiv.sfiDoubleFitCountsLeft{j}  = yhat_left;
            indiv.sfiDoubleFitCountsRight{j} = yhat_right;

            indiv.sfiDoubleFitValid(j) = true;
        end
        fprintf('Scan %d: valid double fits = %d / %d, finite roi1 ratios = %d\n', ...
            i, ...
            nnz(indiv.sfiDoubleFitValid), ...
            Nshots, ...
            nnz(isfinite(indiv.sfiIntegral_roi1_ratio)));        

        % =========================================
        % Plot LEFT/RIGHT amplitudes vs shot index j
        % =========================================
        if isfield(indiv,'sfiIntegralLeft') && isfield(indiv,'sfiIntegralRight')

            Aleft  = NaN(Nshots,1);
            Aright = NaN(Nshots,1);

            for jj = 1:Nshots
                if jj <= numel(indiv.sfiIntegralLeft) && isfinite(indiv.sfiIntegralLeft(jj))
                    Aleft(jj) = indiv.sfiIntegralLeft(jj);
                end
                if jj <= numel(indiv.sfiIntegralRight) && isfinite(indiv.sfiIntegralRight(jj))
                    Aright(jj) = indiv.sfiIntegralRight(jj);
                end
            end

            jvals = (1:Nshots).';
            valid = isfinite(Aleft) | isfinite(Aright);

            if any(valid)
                figure('Name', sprintf('Scan %d: LEFT/RIGHT amplitudes vs shot index', i), 'Color','w');
                hold on; box on;

                plot(jvals(valid), Aleft(valid),  'o-', 'LineWidth', 1.5, 'DisplayName','A_{left}');
                plot(jvals(valid), Aright(valid), 'o-', 'LineWidth', 1.5, 'DisplayName','A_{right}');

                xlabel('Shot index j');
                ylabel('A (counts above baseline)');
                title(sprintf('Scan %d: peak integrals per shot', i));
                legend('Location','best');
                set(gca,'FontSize',14);

                hold off;
            end
        end
        % --- compute average reduced chi^2 for this scan ---
        validChi = isfinite(indiv.sfiDoubleFitChi2Red) & indiv.sfiDoubleFitValid(:).';
        if any(validChi)
            indiv.sfiDoubleFitChi2RedMean = mean(indiv.sfiDoubleFitChi2Red(validChi));
            indiv.sfiDoubleFitChi2RedStd  = std(indiv.sfiDoubleFitChi2Red(validChi));

            fprintf('Scan %d: mean reduced chi^2 = %.4f +/- %.4f over %d valid shots\n', ...
                i, indiv.sfiDoubleFitChi2RedMean, indiv.sfiDoubleFitChi2RedStd, nnz(validChi));
        else
            indiv.sfiDoubleFitChi2RedMean = NaN;
            indiv.sfiDoubleFitChi2RedStd  = NaN;
        end
        % write back scan
        indivDataset{i} = indiv;

        % --- QC PLOTS: plot all fitted shots for this scan ---
        try
            if Nshots > 0

                plotIdx = [];
                for jj = 1:Nshots
                    if jj <= numel(indiv.sfiDoubleFitParams) && ~isempty(indiv.sfiDoubleFitParams{jj})
                        plotIdx(end+1) = jj; %#ok<AGROW>
                    end
                end

                if ~isempty(plotIdx)

                    nValid = numel(plotIdx);
                    nCols = ceil(sqrt(nValid));
                    nRows = ceil(nValid / nCols);

                    figure('Name', sprintf('Scan %d | Double-library QC (%d shots)', i, nValid), ...
                           'Color','w');
                    tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

                    for kk = 1:numel(plotIdx)
                        jj = plotIdx(kk);

                        spec = indiv.mcsSpectra{jj};
                        if isempty(spec) || size(spec,2) < 2
                            continue;
                        end

                        t_plot = spec(:,1);
                        y_plot = spec(:,2);

                        [V_plot, dV_plot] = timeToVoltageAxis(t_plot); %#ok<NASGU>
                        [V_plot, ord] = sort(V_plot);
                        y_plot = y_plot(ord);

                        qfit = indiv.sfiDoubleFitParams{jj};
                        if isempty(qfit)
                            continue;
                        end

                        yhat_left  = qfit(1) .* skewHistModel_A([1, muL, sigmaL, 0, fL, tauL], V_plot);
                        yhat_right = qfit(2) .* skewHistModel_A([1, muR, sigmaR, 0, fR, tauR], V_plot);
                        yhat_tot   = qfit(3) + yhat_left + yhat_right;

                        chi2red = NaN;
                        if isfield(indiv, 'sfiDoubleFitChi2Red') && jj <= numel(indiv.sfiDoubleFitChi2Red)
                            chi2red = indiv.sfiDoubleFitChi2Red(jj);
                        end

                        nexttile;
                        hold on; box on;

                        plot(V_plot, yhat_tot,   'r-',  'LineWidth', 1.5, 'DisplayName','Total fit');
                        plot(V_plot, yhat_left,  '--',  'LineWidth', 1.0, 'DisplayName','Left lib');
                        plot(V_plot, yhat_right, '--',  'LineWidth', 1.0, 'DisplayName','Right lib');
                        plot(V_plot, y_plot,     'ko',  'MarkerSize', 2, 'MarkerFaceColor','k', 'DisplayName','Data');

                        title(sprintf('j = %d | \\chi^2_\\nu = %.3f', jj, chi2red));
                        xlabel('Voltage (V)');
                        ylabel('Counts/bin');
                        set(gca,'FontSize',8);

                        if kk == 1
                            legend('Location','best');
                        end

                        hold off;
                    end

                    if isfield(indiv, 'sfiDoubleFitChi2RedMean') && isfinite(indiv.sfiDoubleFitChi2RedMean)
                        sgtitle(sprintf('Scan %d: all fitted shots | mean \\chi^2_\\nu = %.3f', ...
                            i, indiv.sfiDoubleFitChi2RedMean), 'FontSize', 14);
                    else
                        sgtitle(sprintf('Scan %d: all fitted shots', i), 'FontSize', 14);
                    end
                end
            end
        catch ME
            warning('QC plotting failed for scan %d: %s', i, ME.message);
        end

        % =========================================================
        % Table figure: initial guess vs fitted params for QC shots
        % =========================================================
        % try
        %     nPlot = 5;
        %     if Nshots > 0
        %         plotIdx = unique(round(linspace(1, Nshots, min(nPlot, Nshots))));
        % 
        %         keep = false(size(plotIdx));
        %         for kk = 1:numel(plotIdx)
        %             jj = plotIdx(kk);
        %             keep(kk) = (jj <= numel(indiv.sfiDoubleFitInitParams) && ~isempty(indiv.sfiDoubleFitInitParams{jj})) && ...
        %                        (jj <= numel(indiv.sfiDoubleFitParams)     && ~isempty(indiv.sfiDoubleFitParams{jj}));
        %         end
        %         plotIdx = plotIdx(keep);
        % 
        %         if ~isempty(plotIdx)
        % 
        %             figT = uifigure('Name', sprintf('Scan %d | Guess vs Fit params (QC shots)', i), ...
        %                             'Position',[100 100 720 140*numel(plotIdx)+80]);
        % 
        %             glT = uigridlayout(figT, [numel(plotIdx) 1]);
        %             glT.RowHeight   = repmat({140}, 1, numel(plotIdx));
        %             glT.ColumnWidth = {'1x'};
        %             glT.Padding     = [10 10 10 10];
        %             glT.RowSpacing  = 8;
        % 
        %             paramNames = {'Aleft','Aright','b0'}';
        % 
        %             for kk = 1:numel(plotIdx)
        %                 jj = plotIdx(kk);
        % 
        %                 q0   = indiv.sfiDoubleFitInitParams{jj}(:);
        %                 qfit = indiv.sfiDoubleFitParams{jj}(:);
        % 
        %                 Ttab = table(paramNames, q0, qfit, ...
        %                     'VariableNames', {'Parameter','Guess','Fit'});
        % 
        %                 uit = uitable(glT);
        %                 uit.Data = Ttab;
        %                 uit.ColumnName = Ttab.Properties.VariableNames;
        %                 uit.RowName = {};
        %                 uit.Layout.Row = kk;
        %                 uit.Layout.Column = 1;
        %                 uit.Tooltip = sprintf('Shot j=%d', jj);
        %             end
        %         end
        %     end
        % catch ME
        %     warning('Guess-vs-fit table figure failed for scan %d: %s', i, ME.message);
        % end
    end
end
function [mu_ref, sigma_ref, f_ref, tau_ref] = local_get_template_doubleSkew(libFile, transitionName)
    if ~exist(libFile,'file')
        error('Library file not found: %s', libFile);
    end

    T = readtable(libFile, 'FileType','text', 'Delimiter','\t');

    if ~any(strcmp(T.Properties.VariableNames, 'transition_name'))
        error('Library file missing column: transition_name');
    end

    nameCol = string(T.transition_name);
    target  = string(transitionName);

    idx = find(nameCol == target, 1, 'last');
    if isempty(idx)
        error('Transition "%s" not found in library.', target);
    end

    req = ["mu","sigma","f","tau"];
    for k = 1:numel(req)
        if ~any(strcmp(T.Properties.VariableNames, req(k)))
            error('Library file missing column: %s', req(k));
        end
    end

    mu_ref    = T.mu(idx);
    sigma_ref = abs(T.sigma(idx)) + eps;
    f_ref     = min(max(T.f(idx), 0), 1);
    tau_ref   = abs(T.tau(idx)) + eps;

    if ~isfinite(mu_ref) || ~isfinite(sigma_ref) || ~isfinite(f_ref) || ~isfinite(tau_ref)
        error('Template contains non-finite values for "%s".', target);
    end
end
function export_two_library_summed_by_voltage_csv( ...
    exportFolder, dataSetType, leftTransitionName, rightTransitionName, libFile, ...
    T, X, Xsum, Nsum, ...
    groupedFitInitParams, groupedFitParams, ...
    groupedFitCounts, groupedFitCountsLeft, groupedFitCountsRight, ...
    groupedIntegralLeft, groupedIntegralRight, groupedIntegralTotal, ...
    groupedIntegralLeft_err, groupedIntegralRight_err, groupedB0_err, ...
    groupedChi2, groupedChi2Red, groupedFitValid, ...
    grouped_roi1_ratio, grouped_roi2_ratio, ...
    grouped_roi1_ratio_err, grouped_roi2_ratio_err, ...
    leftTemplate, rightTemplate)

    if ~exist(exportFolder, 'dir')
        mkdir(exportFolder);
    end

    timestamp = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
    filename = timestamp + "_" + dataSetType + ".csv";
    fullpath = fullfile(exportFolder, filename);

    headers = { ...
        'RowType', ...
        'DataSetType', ...
        'Timestamp', ...
        'Section', ...
        'Name', ...
        'Value', ...
        'Units', ...
        'GroupIndex', ...
        'T', ...
        'Nsum', ...
        'PointIndex', ...
        'Time_s', ...
        'Voltage_V', ...
        'Counts_avg_per_bin', ...
        'Counts_sum_per_bin', ...
        'Yerr_avg_per_bin', ...
        'FitTotal_counts_per_bin', ...
        'FitLeft_counts_per_bin', ...
        'FitRight_counts_per_bin', ...
        'Residual_counts_per_bin', ...
        'Aleft', ...
        'Aright', ...
        'Atotal', ...
        'B0', ...
        'Aleft_err', ...
        'Aright_err', ...
        'B0_err', ...
        'Chi2', ...
        'Chi2Red', ...
        'ROI1_ratio', ...
        'ROI2_ratio', ...
        'ROI1_ratio_err', ...
        'ROI2_ratio_err', ...
        'FitValid'};

    C = headers;

    C(end+1,:) = make_two_lib_export_row('metadata', dataSetType, timestamp, 'general', 'leftTransitionName', char(string(leftTransitionName)), '', []);
    C(end+1,:) = make_two_lib_export_row('metadata', dataSetType, timestamp, 'general', 'rightTransitionName', char(string(rightTransitionName)), '', []);
    C(end+1,:) = make_two_lib_export_row('metadata', dataSetType, timestamp, 'general', 'libFile', char(string(libFile)), '', []);
    C(end+1,:) = make_two_lib_export_row('metadata', dataSetType, timestamp, 'general', 'N_groups', numel(T), '', []);

    leftNames = {'muL','sigmaL','fL','tauL'};
    leftUnits = {'V','V','','V'};
    for k = 1:4
        C(end+1,:) = make_two_lib_export_row('library_fixed_parameter', dataSetType, timestamp, 'left_template', leftNames{k}, leftTemplate(k), leftUnits{k}, []);
    end

    rightNames = {'muR','sigmaR','fR','tauR'};
    rightUnits = {'V','V','','V'};
    for k = 1:4
        C(end+1,:) = make_two_lib_export_row('library_fixed_parameter', dataSetType, timestamp, 'right_template', rightNames{k}, rightTemplate(k), rightUnits{k}, []);
    end

    paramNames = {'Aleft','Aright','b0'};
    paramUnits = {'counts','counts','counts_per_bin'};

    for n = 1:numel(T)

        q0 = groupedFitInitParams{n};
        qfit = groupedFitParams{n};

        if ~isempty(q0)
            for k = 1:numel(paramNames)
                C(end+1,:) = make_two_lib_group_param_row( ...
                    'fit_initial_guess', dataSetType, timestamp, ...
                    paramNames{k}, q0(k), paramUnits{k}, ...
                    n, T(n), Nsum(n));
            end
        end

        if ~isempty(qfit)
            for k = 1:numel(paramNames)
                C(end+1,:) = make_two_lib_group_param_row( ...
                    'fit_parameter', dataSetType, timestamp, ...
                    paramNames{k}, qfit(k), paramUnits{k}, ...
                    n, T(n), Nsum(n));
            end
        end

        if isempty(X{n}) || size(X{n},2) < 2
            continue;
        end
        if isempty(groupedFitCounts{n}) || isempty(groupedFitCountsLeft{n}) || isempty(groupedFitCountsRight{n})
            continue;
        end

        specAvg = X{n};
        specSum = Xsum{n};

        t_plot = specAvg(:,1);
        y_avg = specAvg(:,2);
        y_sum = specSum(:,2);

        [V_plot, ~] = timeToVoltageAxis(t_plot);
        [V_plot, ord] = sort(V_plot);

        t_plot = t_plot(ord);
        y_avg = y_avg(ord);
        y_sum = y_sum(ord);

        nContrib = Nsum(n);
        if isfinite(nContrib) && nContrib > 0
            yerr = sqrt(max(y_sum, 1)) ./ nContrib;
        else
            yerr = ones(size(y_avg));
        end

        yhat_tot = groupedFitCounts{n};
        yhat_left = groupedFitCountsLeft{n};
        yhat_right = groupedFitCountsRight{n};
        resid = y_avg - yhat_tot;

        b0fit = NaN;
        if ~isempty(qfit) && numel(qfit) >= 3
            b0fit = qfit(3);
        end

        for k = 1:numel(V_plot)
            C(end+1,:) = { ...
                'data', ...
                char(dataSetType), ...
                char(timestamp), ...
                'summed_by_voltage_grouped_fit', ...
                '', ...
                '', ...
                '', ...
                n, ...
                T(n), ...
                Nsum(n), ...
                k, ...
                t_plot(k), ...
                V_plot(k), ...
                y_avg(k), ...
                y_sum(k), ...
                yerr(k), ...
                yhat_tot(k), ...
                yhat_left(k), ...
                yhat_right(k), ...
                resid(k), ...
                groupedIntegralLeft(n), ...
                groupedIntegralRight(n), ...
                groupedIntegralTotal(n), ...
                b0fit, ...
                groupedIntegralLeft_err(n), ...
                groupedIntegralRight_err(n), ...
                groupedB0_err(n), ...
                groupedChi2(n), ...
                groupedChi2Red(n), ...
                grouped_roi1_ratio(n), ...
                grouped_roi2_ratio(n), ...
                grouped_roi1_ratio_err(n), ...
                grouped_roi2_ratio_err(n), ...
                groupedFitValid(n)};
        end
    end

    writecell(C, fullpath);
    fprintf('Exported double-library summed-by-voltage fit data to:\n%s\n', fullpath);
end


function row = make_two_lib_export_row(rowType, dataSetType, timestamp, section, name, value, units, groupIndex)

    if isempty(groupIndex)
        groupIndex = '';
    end

    row = { ...
        char(rowType), char(dataSetType), char(timestamp), ...
        char(section), char(name), value, char(units), ...
        groupIndex, '', '', '', '', '', '', '', '', '', '', '', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
end


function row = make_two_lib_group_param_row( ...
    rowType, dataSetType, timestamp, paramName, paramValue, units, groupIndex, Tval, NsumVal)

    row = { ...
        char(rowType), char(dataSetType), char(timestamp), ...
        'group_fit_parameter', char(paramName), paramValue, char(units), ...
        groupIndex, Tval, NsumVal, ...
        '', '', '', '', '', '', '', '', '', '', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
end