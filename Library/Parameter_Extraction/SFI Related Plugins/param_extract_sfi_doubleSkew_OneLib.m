function indivDataset = param_extract_sfi_doubleSkew_OneLib(indivDataset, rawSide, libTransitionName, libFile, VminCut)
% Global double-peak mixed fit (RAW + LIB) on counts/bin in voltage domain.
%
% New workflow:
% - Sum all valid mcsSpectra across all indivDataset{i} and all shots j
% - Convert the summed spectrum to voltage
% - Apply a lower-voltage cutoff
% - Fit ONE global summed spectrum
%
% Model:
% - LIB peak: [mu_lib, sigma_lib, f_lib, tau_lib] fixed from library; fit only A_lib
% - RAW peak: fit [A_raw, mu_raw, sigma_raw, f_raw, tau_raw]
% - Shared baseline b0
%
% Fit parameter vector:
%   q = [Araw, mu_raw, sigma_raw, f_raw, tau_raw, Alib, b0]

    if nargin < 4 || isempty(libFile)
        writerPath = which('saveSFIParamLibrary_UI');
        if ~isempty(writerPath)
            libFile = fullfile(fileparts(writerPath), 'sfi_fit_parameter_library.txt');
        else
            here = fileparts(mfilename('fullpath'));
            libFile = fullfile(here, 'sfi_fit_parameter_library.txt');
        end
    end

    if nargin < 5 || isempty(VminCut)
        VminCut = -Inf;
    end
    % ---- Standardized plot-data export controls ----
    savePlotDataBool = true;
    plotDataExportFolder = "E:\SFI_Exports";
    dataSetType = "DoubleFitRawPlusLibrary_GlobalSummed";

    rawSide = lower(strtrim(string(rawSide)));
    if rawSide ~= "left" && rawSide ~= "right"
        error('rawSide must be "left" or "right".');
    end

    [mu_lib, sigma_lib, f_lib, tau_lib] = local_get_template(libFile, libTransitionName);

    [tFull, countsFull, contribList, ok] = sum_mcs_spectra_all_scans(indivDataset);
    if ~ok
        warning('No valid spectra found across all datasets.');
        return;
    end

    [Vfull, dVfull] = timeToVoltageAxis(tFull);
    [Vfull, countsFull, dVfull, tFull] = sort_and_align_vectors(Vfull, countsFull, dVfull, tFull);

    [V, counts, dV, t, keepMask] = apply_voltage_cutoff(Vfull, countsFull, dVfull, tFull, VminCut); %#ok<ASGLU>
    if numel(V) < 10
        warning('Too few points remain after applying voltage cutoff.');
        return;
    end

    plot_raw_global_summed_double_sfi(Vfull, countsFull, VminCut, contribList);

    fitInputs = build_doubleSkew_oneLib_initial_guess(V, counts, dV, rawSide, mu_lib, sigma_lib, f_lib, tau_lib);

    fitStruct = fit_doubleSkew_oneLib_global(V, counts, fitInputs, mu_lib, sigma_lib, f_lib, tau_lib);

    metricStruct = compute_fit_metrics_double_global(V, fitStruct.resid0, fitStruct.residFit, numel(fitStruct.pfit));

    Ttab = make_double_fit_parameter_table(fitInputs.q0, fitStruct.pfit);
    Tchi = make_chi2_table(metricStruct);

    if savePlotDataBool
        export_double_raw_plus_lib_global_csv( ...
            plotDataExportFolder, ...
            dataSetType, ...
            rawSide, ...
            libTransitionName, ...
            libFile, ...
            VminCut, ...
            V, t, dV, counts, ...
            fitInputs, ...
            fitStruct, ...
            metricStruct, ...
            [mu_lib, sigma_lib, f_lib, tau_lib], ...
            size(contribList,1));
    end

    plot_double_fit_qc_window_global( ...
        V, counts, contribList, VminCut, ...
        fitStruct.yhat_tot0, fitStruct.yhat_raw0, fitStruct.yhat_lib0, ...
        fitStruct.yhat_tot, fitStruct.yhat_raw, fitStruct.yhat_lib, ...
        fitStruct.resid0, fitStruct.residFit, ...
        metricStruct, Ttab, Tchi, rawSide, libTransitionName);

    indivDataset = store_global_double_fit_results(indivDataset, ...
        fitStruct, metricStruct, ...
        Vfull, countsFull, tFull, dVfull, ...
        V, counts, t, dV, VminCut, ...
        contribList, rawSide, libTransitionName, ...
        [mu_lib, sigma_lib, f_lib, tau_lib]);
    indivDataset = refit_all_datasets_with_fixed_double_shapes( ...
        indivDataset, VminCut, ...
        fitStruct.mu_raw, fitStruct.sigma_raw, fitStruct.f_raw, fitStruct.tau_raw, ...
        mu_lib, sigma_lib, f_lib, tau_lib);

    saveSFIParamLibrary_UI(indivDataset, []);
end


function [t, counts, contribList, ok] = sum_mcs_spectra_all_scans(indivDataset)

    t = [];
    counts = [];
    contribList = [];
    ok = false;

    for i = 1:numel(indivDataset)
        if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
            continue;
        end
        if ~isfield(indivDataset{i}, 'CounterMCS') || ~isfield(indivDataset{i}, 'mcsSpectra')
            continue;
        end

        nShots = min(indivDataset{i}.CounterMCS, numel(indivDataset{i}.mcsSpectra));

        for j = 1:nShots
            if isempty(indivDataset{i}.mcsSpectra{j})
                continue;
            end

            spec = indivDataset{i}.mcsSpectra{j};
            if size(spec,2) < 2
                continue;
            end

            tj = spec(:,1);
            cj = spec(:,2);

            if isempty(t)
                t = tj(:);
                counts = zeros(size(t));
            else
                if numel(tj) ~= numel(t) || any(abs(tj(:) - t) > 1e-12)
                    warning('Dataset %d shot %d skipped: time axis does not match reference axis.', i, j);
                    continue;
                end
            end

            counts = counts + cj(:);
            contribList(end+1,:) = [i, j]; %#ok<AGROW>
        end
    end

    if isempty(t) || isempty(contribList)
        return;
    end

    ok = true;
end


function [V, counts, dV, t] = sort_and_align_vectors(V, counts, dV, t)
    [V, ord] = sort(V(:));
    counts = counts(ord);
    dV = dV(ord);
    t = t(ord);
end


function [V, counts, dV, t, keepMask] = apply_voltage_cutoff(V, counts, dV, t, VminCut)

    keepMask = (V >= VminCut);

    V = V(keepMask);
    counts = counts(keepMask);
    dV = dV(keepMask);
    t = t(keepMask);
end


function plot_raw_global_summed_double_sfi(Vfull, countsFull, VminCut, contribList)
    figure('Name', 'Global summed double-peak raw spectrum', 'Color', 'w');
    hold on;
    box on;

    plot(Vfull, countsFull, 'k.', 'MarkerSize', 10);

    if isfinite(VminCut)
        xline(VminCut, 'r--', 'LineWidth', 1.2, 'DisplayName', 'V cutoff');
        legend('Data', 'V cutoff', 'Location', 'best');
    end

    xlabel('Voltage (V)');
    ylabel('Counts / bin');
    title(sprintf('Global summed raw V vs Counts (%d total spectra)', size(contribList,1)));
    grid on;
    hold off;
end


function fitInputs = build_doubleSkew_oneLib_initial_guess(V, y, dV, rawSide, mu_lib, sigma_lib, f_lib, tau_lib)

    fitInputs = struct();

    b0_0 = median(y(1:min(20, numel(y))));

    yb = y - b0_0;
    yb(yb < 0) = 0;

    m_lib = skewHistModel_A([1, mu_lib, sigma_lib, 0, f_lib, tau_lib], V);

    denom = (m_lib.' * m_lib);
    if ~isfinite(denom) || denom <= 0
        error('Invalid library model projection denominator.');
    end
    Alib_0 = max(0, (yb.' * m_lib) / denom);

    resid = yb - Alib_0 .* m_lib;
    resid(resid < 0) = 0;

    Vc = V + 0.5 * dV;

    if rawSide == "left"
        mask = (Vc < mu_lib);
    else
        mask = (Vc > mu_lib);
    end

    if any(mask)
        [~, kk] = max(resid(mask));
        Vcand = Vc(mask);
        mu_raw0 = Vcand(kk);
    else
        [~, kk] = max(resid);
        mu_raw0 = Vc(kk);
    end

    if ~isfinite(mu_raw0)
        mu_raw0 = median(Vc);
    end

    sigma_raw0 = estimateSigmaLeftHalfMax(Vc, resid);
    if ~isfinite(sigma_raw0) || sigma_raw0 <= 0
        sigma_raw0 = 0.02 * (max(Vc) - min(Vc));
    end

    f_raw0 = 0.45;
    tau_raw0 = 2.0;
    Araw_0 = max(sum(resid), 0);

    Vmin = min(Vc);
    Vmax = max(Vc);
    Vspan = Vmax - Vmin;
    if Vspan <= 0
        error('Voltage span is not positive.');
    end

    sepV = 0.1;

    if rawSide == "left"
        mu_lb = Vmin;
        mu_ub = min(mu_lib - sepV, Vmax);
    else
        mu_lb = max(mu_lib + sepV, Vmin);
        mu_ub = Vmax;
    end

    if ~(isfinite(mu_lb) && isfinite(mu_ub) && (mu_lb < mu_ub))
        error('Impossible raw-peak ordering constraints. Check rawSide / library peak separation.');
    end

    dVmed = median(dV(isfinite(dV) & dV > 0));
    if ~isfinite(dVmed) || dVmed <= 0
        dVmed = (max(V) - min(V)) / max(numel(V)-1, 1);
    end

    sigMin = max(0.5*dVmed, 1e-6);
    sigMax = 5.0;

    tauMin = 1e-3;
    tauMax = 5.0;

    Atot = sum(yb);
    if ~isfinite(Atot) || Atot <= 0
        Atot = sum(max(y - b0_0, 0));
    end
    Atot = max(Atot, 1);
    Amax = 5 * Atot;

    q0 = [Araw_0, mu_raw0, sigma_raw0, f_raw0, tau_raw0, Alib_0, b0_0];
    lb = [0,      mu_lb,   sigMin,     0.4,    tauMin,   0,      -Inf];
    ub = [Amax,   mu_ub,   sigMax,     1.0,    tauMax,   Amax,   Inf];

    fitInputs.q0 = q0;
    fitInputs.lb = lb;
    fitInputs.ub = ub;
end


function fitStruct = fit_doubleSkew_oneLib_global(V, y, fitInputs, mu_lib, sigma_lib, f_lib, tau_lib)

    q0 = fitInputs.q0;
    lb = fitInputs.lb;
    ub = fitInputs.ub;

    model = @(q, Vaxis) ( ...
          skewHistModel_A([q(1), q(2), q(3), 0, q(4), q(5)], Vaxis) ...
        + skewHistModel_A([q(6), mu_lib, sigma_lib, 0, f_lib, tau_lib], Vaxis) ...
        + q(7) );

    fitStruct = struct();

    fitStruct.yhat_raw0 = skewHistModel_A([q0(1), q0(2), q0(3), 0, q0(4), q0(5)], V);
    fitStruct.yhat_lib0 = skewHistModel_A([q0(6), mu_lib, sigma_lib, 0, f_lib, tau_lib], V);
    fitStruct.yhat_tot0 = fitStruct.yhat_raw0 + fitStruct.yhat_lib0 + q0(7);
    fitStruct.resid0 = fitStruct.yhat_tot0 - y;

    opt = optimoptions('lsqcurvefit', 'Display', 'off');

    [fitStruct.pfit, fitStruct.resnorm, fitStruct.residFit, fitStruct.exitflag, fitStruct.output] = ...
        lsqcurvefit(model, q0, V, y, lb, ub, opt);

    qfit = fitStruct.pfit;

    fitStruct.yhat_raw = skewHistModel_A([qfit(1), qfit(2), qfit(3), 0, qfit(4), qfit(5)], V);
    fitStruct.yhat_lib = skewHistModel_A([qfit(6), mu_lib, sigma_lib, 0, f_lib, tau_lib], V);
    fitStruct.yhat_tot = fitStruct.yhat_raw + fitStruct.yhat_lib + qfit(7);
    fitStruct.residFit = fitStruct.yhat_tot - y;

    fitStruct.Araw = qfit(1);
    fitStruct.mu_raw = qfit(2);
    fitStruct.sigma_raw = qfit(3);
    fitStruct.f_raw = qfit(4);
    fitStruct.tau_raw = qfit(5);
    fitStruct.Alib = qfit(6);
    fitStruct.b0 = qfit(7);

    if fitStruct.mu_raw < mu_lib
        fitStruct.Aleft = fitStruct.Araw;
        fitStruct.Aright = fitStruct.Alib;
    else
        fitStruct.Aleft = fitStruct.Alib;
        fitStruct.Aright = fitStruct.Araw;
    end
end


function metricStruct = compute_fit_metrics_double_global(V, resid0, residFit, p)

    N = numel(V);

    metricStruct = struct();
    metricStruct.N = N;
    metricStruct.p = p;
    metricStruct.chi2_init = sum(resid0.^2) / max(N - p, 1);
    metricStruct.chi2_fit  = sum(residFit.^2) / max(N - p, 1);
end


function Ttab = make_double_fit_parameter_table(q0, qfit)

    paramNames = {'Araw','mu_raw','sigma_raw','f_raw','tau_raw','Alib','b0'}';
    Guess = q0(:);
    Fit = qfit(:);

    Ttab = table(paramNames, Guess, Fit, ...
        'VariableNames', {'Parameter','InitialGuess','FittedValue'});
end


function Tchi = make_chi2_table(metricStruct)

    metricNames = {'chi2_init'; 'chi2_fit'; 'N'; 'p'};
    metricVals  = [metricStruct.chi2_init; metricStruct.chi2_fit; metricStruct.N; metricStruct.p];

    Tchi = table(metricNames, metricVals, ...
        'VariableNames', {'Metric','Value'});
end


function plot_double_fit_qc_window_global( ...
    V, y, contribList, VminCut, ...
    yhat_tot0, yhat_raw0, yhat_lib0, ...
    yhat_tot, yhat_raw, yhat_lib, ...
    resid0, residFit, metricStruct, Ttab, Tchi, rawSide, libTransitionName)

    fig = uifigure('Name', sprintf('Global double-fit QC | %d spectra | chi^2 fit %.3g', ...
        size(contribList,1), metricStruct.chi2_fit), ...
        'Position', [100 100 1000 900]);

    gl = uigridlayout(fig, [4 1]);
    gl.RowHeight   = {'2x','1x','1.2x','0.8x'};
    gl.ColumnWidth = {'1x'};
    gl.Padding     = [10 10 10 10];
    gl.RowSpacing  = 8;

    ax1 = uiaxes(gl);
    hold(ax1, 'on');
    box(ax1, 'on');

    plot(ax1, V, yhat_tot0, '--', 'LineWidth', 1.0);
    plot(ax1, V, yhat_raw0, '--', 'LineWidth', 1.0);
    plot(ax1, V, yhat_lib0, '--', 'LineWidth', 1.0);

    plot(ax1, V, yhat_tot, '-', 'LineWidth', 2.0);
    plot(ax1, V, yhat_raw, '-', 'LineWidth', 1.2);
    plot(ax1, V, yhat_lib, '-', 'LineWidth', 1.2);

    plot(ax1, V, y, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');

    if isfinite(VminCut)
        xline(ax1, VminCut, 'r--', 'LineWidth', 1.0);
    end

    xlabel(ax1, 'Voltage (V)');
    ylabel(ax1, 'Counts / bin');
    title(ax1, sprintf('Global double fit | raw side: %s | lib: %s', rawSide, libTransitionName));
    legend(ax1, {'Init total','Init raw','Init lib','Fit total','Fit raw','Fit lib','Data','V cutoff'}, ...
        'Location', 'best');

    ax2 = uiaxes(gl);
    hold(ax2, 'on');
    box(ax2, 'on');
    plot(ax2, V, resid0, '-', 'LineWidth', 1.0);
    plot(ax2, V, residFit, '-', 'LineWidth', 1.2);
    yline(ax2, 0, 'k--');
    if isfinite(VminCut)
        xline(ax2, VminCut, 'r--', 'LineWidth', 1.0);
    end
    xlabel(ax2, 'Voltage (V)');
    ylabel(ax2, 'Residual (counts/bin)');
    legend(ax2, {'Init resid','Fit resid','Zero','V cutoff'}, 'Location', 'best');

    uit1 = uitable(gl);
    uit1.Data = Ttab;
    uit1.ColumnName = Ttab.Properties.VariableNames;
    uit1.RowName = {};

    uit2 = uitable(gl);
    uit2.Data = Tchi;
    uit2.ColumnName = Tchi.Properties.VariableNames;
    uit2.RowName = {};
end


function indivDataset = store_global_double_fit_results(indivDataset, ...
    fitStruct, metricStruct, ...
    Vfull, countsFull, tFull, dVfull, ...
    Vfit, countsFit, tFit, dVfit, VminCut, ...
    contribList, rawSide, libTransitionName, libTemplate)

    for i = 1:numel(indivDataset)
        if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
            continue;
        end

        indivDataset{i}.globalDoubleFit_mode = 'OneLib';
        indivDataset{i}.globalDoubleFit_rawSide = char(rawSide);
        indivDataset{i}.globalDoubleFit_libTransition = char(libTransitionName);
        indivDataset{i}.globalDoubleFit_libTemplate = libTemplate;

        indivDataset{i}.globalDoubleFitParams = fitStruct.pfit(:).';
        indivDataset{i}.globalDoubleFitCountsTotal = fitStruct.yhat_tot;
        indivDataset{i}.globalDoubleFitCountsRaw = fitStruct.yhat_raw;
        indivDataset{i}.globalDoubleFitCountsLib = fitStruct.yhat_lib;

        indivDataset{i}.globalDoubleFitA_raw = fitStruct.Araw;
        indivDataset{i}.globalDoubleFitA_lib = fitStruct.Alib;
        indivDataset{i}.globalDoubleFitB0 = fitStruct.b0;
        indivDataset{i}.globalDoubleFitA_left = fitStruct.Aleft;
        indivDataset{i}.globalDoubleFitA_right = fitStruct.Aright;

        indivDataset{i}.globalDoubleFitChi2_init = metricStruct.chi2_init;
        indivDataset{i}.globalDoubleFitChi2_fit = metricStruct.chi2_fit;

        indivDataset{i}.globalSummedDoubleSpectraContributors = contribList;

        indivDataset{i}.globalSummedDoubleSpectra_time_full = tFull;
        indivDataset{i}.globalSummedDoubleSpectra_voltage_full = Vfull;
        indivDataset{i}.globalSummedDoubleSpectra_counts_full = countsFull;
        indivDataset{i}.globalSummedDoubleSpectra_dV_full = dVfull;

        indivDataset{i}.globalSummedDoubleSpectra_time_fit = tFit;
        indivDataset{i}.globalSummedDoubleSpectra_voltage_fit = Vfit;
        indivDataset{i}.globalSummedDoubleSpectra_counts_fit = countsFit;
        indivDataset{i}.globalSummedDoubleSpectra_dV_fit = dVfit;

        indivDataset{i}.globalSummedDoubleSpectra_VminCut = VminCut;

        % =================================================
        % Compatibility fields for saveSFIParamLibrary_UI
        % Save RAW component as a standard 6-parameter peak
        % [A, mu, sigma, b0, f, tau]
        % =================================================
        indivDataset{i}.globalSkewFitParams_summed = [ ...
            fitStruct.Araw, ...
            fitStruct.mu_raw, ...
            fitStruct.sigma_raw, ...
            fitStruct.b0, ...
            fitStruct.f_raw, ...
            fitStruct.tau_raw ...
        ];

        indivDataset{i}.globalSummedSFISpectra_dV = median(dVfit);
        indivDataset{i}.globalSummedSFISpectraContributors = contribList;
    end
end


function [mu_ref, sigma_ref, f_ref, tau_ref] = local_get_template(libFile, transitionName)

    if ~exist(libFile,'file')
        error('Library file not found: %s', libFile);
    end

    T = readtable(libFile, 'FileType', 'text', 'Delimiter', '\t');

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
function indivDataset = refit_all_datasets_with_fixed_double_shapes( ...
    indivDataset, VminCut, ...
    mu_raw, sigma_raw, f_raw, tau_raw, ...
    mu_lib, sigma_lib, f_lib, tau_lib)

    for i = 1:numel(indivDataset)
        if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
            continue;
        end
        if ~isfield(indivDataset{i}, 'CounterMCS') || ~isfield(indivDataset{i}, 'mcsSpectra')
            continue;
        end

        indiv = indivDataset{i};
        nShots = min(indiv.CounterMCS, numel(indiv.mcsSpectra));

        indiv.fixedDoubleRefitParams = cell(1, nShots);       % [Araw, Alib, b0]
        indiv.fixedDoubleRefitChi2Red = nan(1, nShots);
        indiv.fixedDoubleRefitCountsTotal = cell(1, nShots);
        indiv.fixedDoubleRefitCountsRaw = cell(1, nShots);
        indiv.fixedDoubleRefitCountsLib = cell(1, nShots);
        indiv.fixedDoubleRefitValid = false(1, nShots);

        plotIdx = [];

        for j = 1:nShots
            if isempty(indiv.mcsSpectra{j})
                continue;
            end

            spec = indiv.mcsSpectra{j};
            if size(spec,2) < 2
                continue;
            end

            t = spec(:,1);
            y = spec(:,2);

            if numel(t) < 5
                continue;
            end

            t = t(:);
            y = y(:);

            [V, dV] = timeToVoltageAxis(t);
            [V, y, dV, t] = sort_and_align_vectors(V, y, dV, t);
            [V, y, dV, t] = apply_voltage_cutoff(V, y, dV, t, VminCut); %#ok<ASGLU>

            if numel(V) < 5 || any(~isfinite(V)) || any(~isfinite(y))
                continue;
            end

            refitStruct = fit_single_spectrum_fixed_double_shapes( ...
                V, y, mu_raw, sigma_raw, f_raw, tau_raw, ...
                mu_lib, sigma_lib, f_lib, tau_lib);

            indiv.fixedDoubleRefitParams{j} = refitStruct.qfit(:).';
            indiv.fixedDoubleRefitChi2Red(j) = refitStruct.chi2red;
            indiv.fixedDoubleRefitCountsTotal{j} = refitStruct.yhat_tot;
            indiv.fixedDoubleRefitCountsRaw{j} = refitStruct.yhat_raw;
            indiv.fixedDoubleRefitCountsLib{j} = refitStruct.yhat_lib;
            indiv.fixedDoubleRefitValid(j) = true;

            plotIdx(end+1) = j; %#ok<AGROW>
        end

        if ~isempty(plotIdx)
            make_fixed_double_refit_dataset_figure(indiv, i, plotIdx, VminCut);
        end

        validChi = isfinite(indiv.fixedDoubleRefitChi2Red);
        if any(validChi)
            indiv.fixedDoubleRefitChi2RedMean = mean(indiv.fixedDoubleRefitChi2Red(validChi));
        else
            indiv.fixedDoubleRefitChi2RedMean = NaN;
        end

        indivDataset{i} = indiv;
    end
end


function refitStruct = fit_single_spectrum_fixed_double_shapes( ...
    V, y, mu_raw, sigma_raw, f_raw, tau_raw, ...
    mu_lib, sigma_lib, f_lib, tau_lib)

    % Fixed shapes, only fit [Araw, Alib, b0]
    m_raw = skewHistModel_A([1, mu_raw, sigma_raw, 0, f_raw, tau_raw], V);
    m_lib = skewHistModel_A([1, mu_lib, sigma_lib, 0, f_lib, tau_lib], V);

    b0_0 = median(y(1:min(20, numel(y))));
    X = [m_raw(:), m_lib(:), ones(numel(V),1)];
    q0 = X \ y(:);

    if numel(q0) < 3 || any(~isfinite(q0))
        q0 = [0; 0; b0_0];
    end

    q0(1) = max(q0(1), 0);
    q0(2) = max(q0(2), 0);
    q0(3) = b0_0;

    yb = y - b0_0;
    yb(yb < 0) = 0;
    Atot = sum(yb);
    if ~isfinite(Atot) || Atot <= 0
        Atot = sum(max(y - b0_0, 0));
    end
    Atot = max(Atot, 1);
    Amax = 5 * Atot;

    lb = [0,    0,    -Inf];
    ub = [Amax, Amax,  Inf];

    model = @(q, Vaxis) ( ...
        q(1) .* skewHistModel_A([1, mu_raw, sigma_raw, 0, f_raw, tau_raw], Vaxis) + ...
        q(2) .* skewHistModel_A([1, mu_lib, sigma_lib, 0, f_lib, tau_lib], Vaxis) + ...
        q(3));

    opt = optimoptions('lsqcurvefit', 'Display', 'off');
    qfit = lsqcurvefit(model, q0(:).', V, y, lb, ub, opt);

    yhat_raw = qfit(1) .* skewHistModel_A([1, mu_raw, sigma_raw, 0, f_raw, tau_raw], V);
    yhat_lib = qfit(2) .* skewHistModel_A([1, mu_lib, sigma_lib, 0, f_lib, tau_lib], V);
    yhat_tot = yhat_raw + yhat_lib + qfit(3);

    resid = yhat_tot - y;
    N = numel(V);
    p = 3;
    chi2red = sum(resid.^2) / max(N - p, 1);

    refitStruct = struct();
    refitStruct.qfit = qfit;
    refitStruct.yhat_raw = yhat_raw;
    refitStruct.yhat_lib = yhat_lib;
    refitStruct.yhat_tot = yhat_tot;
    refitStruct.resid = resid;
    refitStruct.chi2red = chi2red;
end


function make_fixed_double_refit_dataset_figure(indiv, scanIdx, plotIdx, VminCut)

    nValid = numel(plotIdx);
    nCols = ceil(sqrt(nValid));
    nRows = ceil(nValid / nCols);

    validChi = isfinite(indiv.fixedDoubleRefitChi2Red);
    if any(validChi)
        meanChi = mean(indiv.fixedDoubleRefitChi2Red(validChi));
    else
        meanChi = NaN;
    end

    figure('Name', sprintf('Scan %d | fixed-shape double refits', scanIdx), ...
           'Color', 'w');
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for kk = 1:nValid
        j = plotIdx(kk);

        spec = indiv.mcsSpectra{j};
        t = spec(:,1);
        y = spec(:,2);

        [V, dV] = timeToVoltageAxis(t);
        [V, y, dV, t] = sort_and_align_vectors(V, y, dV, t);
        [V, y, dV, t] = apply_voltage_cutoff(V, y, dV, t, VminCut); %#ok<ASGLU>

        yhat_tot = indiv.fixedDoubleRefitCountsTotal{j};
        yhat_raw = indiv.fixedDoubleRefitCountsRaw{j};
        yhat_lib = indiv.fixedDoubleRefitCountsLib{j};
        chi2red = indiv.fixedDoubleRefitChi2Red(j);

        nexttile;
        hold on;
        box on;

        plot(V, yhat_tot, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Total fit');
        plot(V, yhat_raw, '--', 'LineWidth', 1.0, 'DisplayName', 'Raw');
        plot(V, yhat_lib, '--', 'LineWidth', 1.0, 'DisplayName', 'Lib');
        plot(V, y, 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Data');

        if isfinite(VminCut)
            xline(VminCut, 'r--', 'LineWidth', 1.0);
        end

        title(sprintf('j = %d | \\chi^2_\\nu = %.3f', j, chi2red));
        xlabel('Voltage (V)');
        ylabel('Counts/bin');
        set(gca, 'FontSize', 8);

        if kk == 1
            legend('Location', 'best');
        end

        hold off;
    end

    sgtitle(sprintf('Scan %d | mean \\chi^2_\\nu = %.3f', scanIdx, meanChi), ...
        'FontSize', 14);
end
function export_double_raw_plus_lib_global_csv( ...
    exportFolder, dataSetType, rawSide, libTransitionName, libFile, VminCut, ...
    V, t, dV, counts, fitInputs, fitStruct, metricStruct, libTemplate, nContrib)

    if ~exist(exportFolder, 'dir')
        mkdir(exportFolder);
    end

    timestamp = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
    filename = timestamp + "_" + dataSetType + ".csv";
    fullpath = fullfile(exportFolder, filename);

    mu_lib    = libTemplate(1);
    sigma_lib = libTemplate(2);
    f_lib     = libTemplate(3);
    tau_lib   = libTemplate(4);

    headers = { ...
        'RowType', ...
        'DataSetType', ...
        'Timestamp', ...
        'Section', ...
        'Name', ...
        'Value', ...
        'Units', ...
        'PointIndex', ...
        'Time_s', ...
        'Voltage_V', ...
        'dV_V', ...
        'Counts_per_bin', ...
        'InitTotal_counts_per_bin', ...
        'InitRaw_counts_per_bin', ...
        'InitLib_counts_per_bin', ...
        'FitTotal_counts_per_bin', ...
        'FitRaw_counts_per_bin', ...
        'FitLib_counts_per_bin', ...
        'InitResidual_counts_per_bin', ...
        'FitResidual_counts_per_bin'};

    C = headers;

    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'rawSide', char(rawSide), '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'libTransitionName', char(string(libTransitionName)), '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'libFile', char(string(libFile)), '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'VminCut', VminCut, 'V', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'N_contributing_spectra', nContrib, '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'N_points', metricStruct.N, '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'general', 'p', metricStruct.p, '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'fit_metric', 'chi2_init', metricStruct.chi2_init, '', []);
    C(end+1,:) = make_double_export_row('metadata', dataSetType, timestamp, 'fit_metric', 'chi2_fit', metricStruct.chi2_fit, '', []);

    paramNames = {'Araw','mu_raw','sigma_raw','f_raw','tau_raw','Alib','b0'};
    units      = {'counts','V','V','','V','counts','counts_per_bin'};

    for k = 1:numel(paramNames)
        C(end+1,:) = make_double_export_row('fit_initial_guess', dataSetType, timestamp, 'fit_parameter', paramNames{k}, fitInputs.q0(k), units{k}, []);
        C(end+1,:) = make_double_export_row('fit_parameter', dataSetType, timestamp, 'fit_parameter', paramNames{k}, fitStruct.pfit(k), units{k}, []);
    end

    C(end+1,:) = make_double_export_row('library_fixed_parameter', dataSetType, timestamp, 'lib_template', 'mu_lib', mu_lib, 'V', []);
    C(end+1,:) = make_double_export_row('library_fixed_parameter', dataSetType, timestamp, 'lib_template', 'sigma_lib', sigma_lib, 'V', []);
    C(end+1,:) = make_double_export_row('library_fixed_parameter', dataSetType, timestamp, 'lib_template', 'f_lib', f_lib, '', []);
    C(end+1,:) = make_double_export_row('library_fixed_parameter', dataSetType, timestamp, 'lib_template', 'tau_lib', tau_lib, 'V', []);

    for k = 1:numel(V)
        C(end+1,:) = { ...
            'data', char(dataSetType), char(timestamp), ...
            'global_summed_double_fit', '', '', '', ...
            k, t(k), V(k), dV(k), counts(k), ...
            fitStruct.yhat_tot0(k), ...
            fitStruct.yhat_raw0(k), ...
            fitStruct.yhat_lib0(k), ...
            fitStruct.yhat_tot(k), ...
            fitStruct.yhat_raw(k), ...
            fitStruct.yhat_lib(k), ...
            fitStruct.resid0(k), ...
            fitStruct.residFit(k)};
    end

    writecell(C, fullpath);
    fprintf('Exported global double raw+library fit data to:\n%s\n', fullpath);
end


function row = make_double_export_row(rowType, dataSetType, timestamp, section, name, value, units, pointIndex)

    if isempty(pointIndex)
        pointIndex = '';
    end

    row = { ...
        char(rowType), char(dataSetType), char(timestamp), ...
        char(section), char(name), value, char(units), ...
        pointIndex, '', '', '', '', '', '', '', '', '', '', '', ''};
end