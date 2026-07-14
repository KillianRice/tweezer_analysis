function indivDataset = param_extract_sfi_integral_fitting(analyVar, indivDataset, avgDataset)

    disp('Entered param_extract_sfi_integral_fitting');
    if nargin < 3
        avgDataset = [];
    end

    [t, counts, contribList, ok] = sum_mcs_spectra_all_scans(indivDataset, analyVar);
    if ~ok
        warning('No valid spectra found across all datasets.');
        return;
    end
    % --- standardized plot-data export controls ---
    savePlotDataBool = true;
    plotDataExportFolder = "E:\SFI_Exports";
    dataSetType = " SingleFitRaw_GlobalSummedSFI";


    [V, dVvec] = timeToVoltageAxis(t);
    dV = median(dVvec);

    % --- apply lower voltage cutoff ---
    VminCut = 30;   % example value, change as needed

    keepMask = V >= VminCut;

    V = V(keepMask);
    counts = counts(keepMask);
    t = t(keepMask);

    if numel(V) < 10
        warning('Too few points remain after voltage cutoff.');
        return;
    end

    plot_raw_global_summed_sfi(V, counts, contribList);

    model = @(p, xx) skewHistModel_A(p, xx);

    [p0, lb, ub] = build_initial_guess_and_bounds(V, counts);

    fitStruct = fit_summed_sfi_spectrum(model, p0, lb, ub, V, counts);

metricStruct = compute_fit_metrics(V, fitStruct.resid0, fitStruct.residFit);

Ttab = make_fit_parameter_table(p0, fitStruct.pfit);
Tchi = make_chi2_table(metricStruct);

% ============================================
% export plot / fit data
% ============================================
if savePlotDataBool
    export_sfi_single_fit_csv( ...
        plotDataExportFolder, ...
        dataSetType, ...
        V, t, counts, ...
        fitStruct, ...
        metricStruct, ...
        Ttab, Tchi, ...
        contribList, ...
        dV, ...
        VminCut);
end

plot_fit_qc_window_global(V, counts, contribList, ...
    fitStruct.rhat0, fitStruct.rhatFit, ...
    fitStruct.resid0, fitStruct.residFit, ...
    metricStruct, Ttab, Tchi);

    indivDataset = store_global_summed_fit_results(indivDataset, dV, fitStruct.pfit, counts, t, V, contribList);

    saveSFIParamLibrary_UI(indivDataset, analyVar);
end


function [t, counts, contribList, ok] = sum_mcs_spectra_all_scans(indivDataset, analyVar)

    t = [];
    counts = [];
    contribList = [];
    ok = false;

    nDatasets = min(analyVar.numBasenamesAtom, numel(indivDataset));

    for i = 1:nDatasets
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
                if numel(tj) ~= numel(t) || any(tj(:) ~= t)
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


function plot_raw_global_summed_sfi(V, counts, contribList)
    figure;
    plot(V, counts, 'k.', 'MarkerSize', 10);
    xlabel('Voltage (V)');
    ylabel('Counts / bin');
    title(sprintf('Global summed raw V vs Counts (%d total spectra)', size(contribList,1)));
    grid on;
    box on;
end


function [p0, lb, ub] = build_initial_guess_and_bounds(V, counts)

    b_0 = median(counts(1:min(20, numel(counts))));

    A_0 = counts - b_0;
    A_0(A_0 < 0) = 0;
    A_0 = sum(A_0);

    [~, ipk] = max(counts);
    x_0 = V(ipk);

    sigma_0 = estimateSigmaLeftHalfMax(V, max(counts - b_0, 0));
    f_0 = 0.45;
    tau_0 = sigma_0 * 7;

    p0 = [A_0, x_0, sigma_0, b_0, f_0, tau_0];

    lb = [0, min(V), 0, 0, 0.3, 0];
    ub = [Inf, max(V), Inf, Inf, 1, Inf];
end


function fitStruct = fit_summed_sfi_spectrum(model, p0, lb, ub, V, counts)

    fitStruct = struct();

    fitStruct.rhat0 = model(p0, V);
    fitStruct.resid0 = fitStruct.rhat0 - counts;

    [fitStruct.pfit, fitStruct.resnorm, fitStruct.residFit, fitStruct.exitflag, fitStruct.output] = ...
        lsqcurvefit(model, p0, V, counts, lb, ub);

    fitStruct.rhatFit = model(fitStruct.pfit, V);
    fitStruct.residFit = fitStruct.rhatFit - counts;

    fprintf('max(counts)=%.3g, max(rhat0)=%.3g, max(rhatFit)=%.3g\n', ...
        max(counts), max(fitStruct.rhat0), max(fitStruct.rhatFit));
end


function metricStruct = compute_fit_metrics(V, resid0, residFit)

    N = numel(V);
    p = 6;

    metricStruct = struct();
    metricStruct.N = N;

    metricStruct.chi2_init = sum(resid0.^2) / (N - p);
    metricStruct.chi2_fit  = sum(residFit.^2) / (N - p);
end


function Ttab = make_fit_parameter_table(p0, pfit)
    paramNames   = {'A','mu','sigma','b0','f','tau'}';
    InitialGuess = p0(:);
    FittedValue  = pfit(:);

    Ttab = table(paramNames, InitialGuess, FittedValue, ...
        'VariableNames', {'Parameter','InitialGuess','FittedValue'});
end


function Tchi = make_chi2_table(metricStruct)
    chiNames = {'chi2_init'; 'chi2_fit'};
    chiVals  = [metricStruct.chi2_init; metricStruct.chi2_fit];

    Tchi = table(chiNames, chiVals, ...
        'VariableNames', {'Metric','Value'});
end


function plot_fit_qc_window_global(V, counts, contribList, rhat0, rhatFit, resid0, residFit, metricStruct, Ttab, Tchi)

    fig = uifigure('Name', sprintf('Global summed fit | %d spectra | chi^2 init %.3g | chi^2 fit %.3g', ...
        size(contribList,1), metricStruct.chi2_init, metricStruct.chi2_fit), ...
        'Position', [100 100 950 820]);

    gl = uigridlayout(fig, [4 1]);
    gl.RowHeight   = {'2x','1x','1.2x','0.6x'};
    gl.ColumnWidth = {'1x'};
    gl.Padding     = [10 10 10 10];
    gl.RowSpacing  = 8;

    ax1 = uiaxes(gl);
    hold(ax1, 'on');
    box(ax1, 'on');
    plot(ax1, V, rhat0, '-', 'LineWidth', 1.2);
    plot(ax1, V, rhatFit, '-', 'LineWidth', 2.0);
    plot(ax1, V, counts, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    xlabel(ax1, 'Voltage (V)');
    ylabel(ax1, 'Counts / bin');
    legend(ax1, {'Initial','Fit','Data'}, 'Location', 'best');

    ax2 = uiaxes(gl);
    hold(ax2, 'on');
    box(ax2, 'on');
    plot(ax2, V, resid0, '-', 'LineWidth', 1.0);
    plot(ax2, V, residFit, '-', 'LineWidth', 1.2);
    yline(ax2, 0, 'k--');
    xlabel(ax2, 'Voltage (V)');
    ylabel(ax2, 'Residual (counts/bin)');
    legend(ax2, {'Init resid','Fit resid'}, 'Location', 'best');

    uit = uitable(gl);
    uit.Data = Ttab;
    uit.ColumnName = Ttab.Properties.VariableNames;
    uit.RowName = {};

    uit2 = uitable(gl);
    uit2.Data = Tchi;
    uit2.ColumnName = Tchi.Properties.VariableNames;
    uit2.RowName = {};
end


function indivDataset = store_global_summed_fit_results(indivDataset, dV, pfit, counts, t, V, contribList)

    for i = 1:numel(indivDataset)
        if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
            continue;
        end

        indivDataset{i}.globalSummedSFISpectra_dV = dV;
        indivDataset{i}.globalSkewFitParams_summed = pfit;
        indivDataset{i}.globalSummedSFISpectra_counts = counts;
        indivDataset{i}.globalSummedSFISpectra_time = t;
        indivDataset{i}.globalSummedSFISpectra_voltage = V;
        indivDataset{i}.globalSummedSFISpectraContributors = contribList;
    end
end
function export_sfi_single_fit_csv(exportFolder, dataSetType, ...
    V, t, counts, fitStruct, metricStruct, Ttab, Tchi, ...
    contribList, dV, VminCut)

    if ~exist(exportFolder, 'dir')
        mkdir(exportFolder);
    end

    timestamp = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
    filename = timestamp + "_" + dataSetType + ".csv";
    fullpath = fullfile(exportFolder, filename);

    headers = {'RowType','DataSetType','Timestamp','Section','Name','Value','Units', ...
               'Index','Voltage_V','Time_s','Counts_per_bin', ...
               'InitialFit_counts_per_bin','FittedFit_counts_per_bin', ...
               'InitialResidual_counts_per_bin','FittedResidual_counts_per_bin', ...
               'DatasetIndex','ShotIndex'};

    C = headers;

    % -------------------------
    % Metadata rows
    % -------------------------
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"general","model","skewHistModel_A","",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"general","dV",dV,"V",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"general","VminCut",VminCut,"V",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"general","N_points",metricStruct.N,"",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"general","N_contributing_spectra",size(contribList,1),"",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"fit_metric","chi2_init",metricStruct.chi2_init,"",[]);
    C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"fit_metric","chi2_fit",metricStruct.chi2_fit,"",[]);

    for k = 1:height(Ttab)
        pname = string(Ttab.Parameter{k});
        C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"initial_guess",pname,Ttab.InitialGuess(k),"",[]);
        C(end+1,:) = make_csv_row("metadata",dataSetType,timestamp,"fitted_parameter",pname,Ttab.FittedValue(k),"",[]);
    end

    % -------------------------
    % Main plot/reconstruction rows
    % -------------------------
    for k = 1:numel(V)
        row = { ...
            'data', char(dataSetType), char(timestamp), 'global_summed_sfi', '', '', '', ...
            k, V(k), t(k), counts(k), ...
            fitStruct.rhat0(k), fitStruct.rhatFit(k), ...
            fitStruct.resid0(k), fitStruct.residFit(k), ...
            '', ''};
        C(end+1,:) = row;
    end

    writecell(C, fullpath);

    fprintf('Exported plot/fit data to:\n%s\n', fullpath);
end


function row = make_csv_row(rowType,dataSetType,timestamp,section,name,value,units,idx)
    if isempty(idx)
        idx = '';
    end

    row = { ...
        char(rowType), char(dataSetType), char(timestamp), char(section), char(name), value, char(units), ...
        idx, '', '', '', '', '', '', '', '', ''};
end