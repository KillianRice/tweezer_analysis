function indivDataset = sfi_fit_allshots_library_singlepeak(indivDataset, transitionName, libFile)
% Fit all mcsSpectra in ALL indivDataset{i} entries using a single-peak library template:
%   - mu, sigma, f, tau fixed from library (shape parameters, in mV)
%   - fit ONLY (A, b0) per shot by linear least squares
%   - fit is performed on raw histogram counts (counts/bin) vs voltage V (mV)
%
% Model for expected counts/bin:
%   chat(V) = b0 + A * m(V; mu,sigma,f,tau)
%
% OUTPUTS added/overwritten in each indivDataset{i}:
%   .sfiIntegral(j)     = Afit  numeric vector
%   .sfiFitA(j)         = Afit
%   .sfiFitB0(j)        = b0fit
%   .sfiFitParams{j}    = [Afit mu sigma b0fit f tau]
%   .sfiFitTransition
%   .sfiFitTemplate     = [mu sigma f tau]

    % ---- Library path ----
    if nargin < 3 || isempty(libFile)
        here = fileparts(mfilename('fullpath'));
        libFile = fullfile(here, 'sfi_fit_parameter_library.txt');
    end

    % ---- Standardized plot-data export controls ----
    savePlotDataBool = false;
    plotDataExportFolder = "E:\SFI_Exports";
    dataSetType = "SingleFitLibrary_AllShots";

    % ---- Load template once ----
    [mu_ref, sigma_ref, f_ref, tau_ref] = local_get_template(libFile, transitionName);

    % ---- Initialize export container ----
    exportRows = {};
    exportTimestamp = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));

    % ---- Loop over all scans i ----
    for i = 1:numel(indivDataset)

        indiv = indivDataset{i};

        if isempty(indiv) || ~isstruct(indiv)
            continue;
        end
        if ~isfield(indiv,'CounterMCS') || indiv.CounterMCS < 1
            continue;
        end
        if ~isfield(indiv,'mcsSpectra') || isempty(indiv.mcsSpectra)
            continue;
        end

        % Store template metadata
        indiv.sfiFitTransition = char(transitionName);
        indiv.sfiFitTemplate   = [mu_ref, sigma_ref, f_ref, tau_ref];

        % ------------------------------------------------------------
        % IMPORTANT:
        % sfiIntegral must be numeric, not cell, so getxy / average_plot
        % can use ydata_clean{basename}(j) directly.
        % ------------------------------------------------------------
        indiv.sfiIntegral  = nan(1, indiv.CounterMCS);
        indiv.sfiFitA      = nan(1, indiv.CounterMCS);
        indiv.sfiFitB0     = nan(1, indiv.CounterMCS);
        indiv.sfiFitParams = cell(1, indiv.CounterMCS);

        % ---- Fit all shots j ----
        for j = 1:indiv.CounterMCS

            if j > numel(indiv.mcsSpectra) || isempty(indiv.mcsSpectra{j})
                continue;
            end

            spec = indiv.mcsSpectra{j};

            if size(spec,2) < 2
                continue;
            end

            t = spec(:,1);
            y = spec(:,2);   % raw counts/bin

            if numel(t) < 3 || numel(y) ~= numel(t)
                continue;
            end
            if any(~isfinite(t)) || any(~isfinite(y))
                continue;
            end

            % --- Convert time axis -> voltage axis (mV) ---
            [V, dV] = timeToVoltageAxis(t);   % V in mV, dV in mV

            % --- Build shape vector m(V) for A=1, b0=0 using fixed template ---
            m = local_shape_counts_per_bin_V(V, dV, mu_ref, sigma_ref, f_ref, tau_ref);

            % --- Linear least squares for [b0; A] ---
            % y ≈ b0*1 + A*m
            X = [ones(size(V)), m];
            x = X \ y;

            b0fit = x(1);
            Afit  = x(2);

            % Simple sanity clamps
            if ~isfinite(b0fit)
                b0fit = 0;
            end
            if ~isfinite(Afit)
                Afit = 0;
            end
            if Afit < 0
                Afit = 0;
            end

            % Predicted counts/bin on the same bins
            yhat = b0fit + Afit .* m;

            % ---- Add this shot to export data ----
            if savePlotDataBool
                residual = yhat - y;

                exportRows = append_library_singlefit_export_rows( ...
                    exportRows, ...
                    dataSetType, ...
                    exportTimestamp, ...
                    transitionName, ...
                    libFile, ...
                    i, j, ...
                    t, V, dV, y, m, yhat, residual, ...
                    Afit, b0fit, ...
                    mu_ref, sigma_ref, f_ref, tau_ref);

                exportRows = append_library_singlefit_param_rows( ...
                    exportRows, ...
                    dataSetType, ...
                    exportTimestamp, ...
                    transitionName, ...
                    libFile, ...
                    i, j, ...
                    Afit, mu_ref, sigma_ref, b0fit, f_ref, tau_ref);
            end

            % Store results
            indiv.sfiIntegral(j)  = Afit;
            indiv.sfiFitA(j)      = Afit;
            indiv.sfiFitB0(j)     = b0fit;
            indiv.sfiFitParams{j} = [Afit, mu_ref, sigma_ref, b0fit, f_ref, tau_ref];

            % Optional: store trace if you want later
            % indiv.sfiFitYhat{j} = yhat;
        end

        % Write back
        indivDataset{i} = indiv;

        % ---- Plot Afit vs shot index j ----
        if any(isfinite(indiv.sfiFitA))
            figure('Name',sprintf('Scan %d: Afit vs shot index', i));
            plot(1:numel(indiv.sfiFitA), indiv.sfiFitA, 'ko-', ...
                'LineWidth', 1.5, 'MarkerSize', 6);
            xlabel('Shot index j');
            ylabel('A (total counts above baseline)');
            title(sprintf('Scan %d: Library-fit A vs shot index', i));
            box on;
            grid on;
        end

        % ---- Plot all valid shot fits in one tiled window ----
        validShots = find(isfinite(indiv.sfiFitA));
        nValid = numel(validShots);

        if nValid > 0
            nCols = ceil(sqrt(nValid));
            nRows = ceil(nValid / nCols);

            figure('Name', sprintf('Scan %d - All Shot Fits (Library)', i));
            tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

            for idx = 1:nValid
                j = validShots(idx);

                spec = indiv.mcsSpectra{j};
                t = spec(:,1);
                counts = spec(:,2);

                [V, dV] = timeToVoltageAxis(t);

                pfit = indiv.sfiFitParams{j};  % [A mu sigma b0 f tau]
                Afit  = pfit(1);
                mu    = pfit(2);
                sigma = pfit(3);
                b0fit = pfit(4);
                f     = pfit(5);
                tau   = pfit(6);

                m = local_shape_counts_per_bin_V(V, dV, mu, sigma, f, tau);
                fitCounts = b0fit + Afit .* m;

                nexttile;
                hold on;
                box on;
                plot(V, counts, 'k.', 'MarkerSize', 4);
                plot(V, fitCounts, 'r-', 'LineWidth', 1.2);
                title(sprintf('j=%d', j));
                set(gca,'FontSize',8);
                hold off;
            end

            sgtitle(sprintf('Scan %d: All Shot Fits (Library template = %s)', ...
                i, string(transitionName)), 'FontSize', 14);
        end
    end

    if savePlotDataBool
        export_library_singlefit_csv( ...
            plotDataExportFolder, ...
            dataSetType, ...
            exportTimestamp, ...
            transitionName, ...
            libFile, ...
            mu_ref, sigma_ref, f_ref, tau_ref, ...
            exportRows);
    end
end

% ===================== helpers =====================

function [mu_ref, sigma_ref, f_ref, tau_ref] = local_get_template(libFile, transitionName)

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
    sigma_ref = T.sigma(idx);
    f_ref     = T.f(idx);
    tau_ref   = T.tau(idx);

    if ~isfinite(mu_ref) || ~isfinite(sigma_ref) || ...
            ~isfinite(f_ref) || ~isfinite(tau_ref)
        error('Template contains non-finite values for "%s".', target);
    end

    sigma_ref = abs(sigma_ref) + eps;
    tau_ref   = abs(tau_ref)   + eps;
    f_ref     = min(max(f_ref, 0), 1);
end

function m = local_shape_counts_per_bin_V(V, dV, mu, sigma, f, tau)
% Build m_i = p(Vc_i) * dV_i using bin centers and per-bin widths.
% V, dV in mV. mu, sigma, tau in mV.

    V  = V(:);
    dV = dV(:);

    sigma = abs(sigma) + eps;
    tau   = abs(tau)   + eps;
    f     = min(max(f,0),1);

    % Bin centers
    Vc = V + 0.5*dV;

    % Normalized Gaussian pdf, area = 1
    g = (1./(sigma*sqrt(2*pi))) .* exp(-0.5 .* ((Vc - mu)./sigma).^2);

    % Normalized one-sided exponential pdf, area = 1
    u = Vc - mu;
    e = (1./tau) .* exp(-max(u,0)./tau) .* double(u >= 0);

    pdf = (1 - f).*g + f.*e;

    % Counts/bin contribution for A = 1
    m = pdf .* dV;
end

function exportRows = append_library_singlefit_export_rows( ...
    exportRows, dataSetType, exportTimestamp, transitionName, libFile, ...
    scanIndex, shotIndex, t, V, dV, counts, shape_m, fitCounts, residual, ...
    Afit, b0fit, mu_ref, sigma_ref, f_ref, tau_ref)

    for k = 1:numel(V)
        exportRows(end+1,:) = { ...
            'data', ...
            char(dataSetType), ...
            char(exportTimestamp), ...
            char(string(transitionName)), ...
            char(string(libFile)), ...
            scanIndex, ...
            shotIndex, ...
            k, ...
            t(k), ...
            V(k), ...
            dV(k), ...
            counts(k), ...
            shape_m(k), ...
            fitCounts(k), ...
            residual(k), ...
            Afit, ...
            b0fit, ...
            mu_ref, ...
            sigma_ref, ...
            f_ref, ...
            tau_ref};
    end
end

function export_library_singlefit_csv( ...
    exportFolder, dataSetType, exportTimestamp, transitionName, libFile, ...
    mu_ref, sigma_ref, f_ref, tau_ref, exportRows)

    if ~exist(exportFolder, 'dir')
        mkdir(exportFolder);
    end

    filename = exportTimestamp + "_" + dataSetType + ".csv";
    fullpath = fullfile(exportFolder, filename);

    headers = { ...
        'RowType', ...
        'DataSetType', ...
        'Timestamp', ...
        'TransitionName', ...
        'LibraryFile', ...
        'ScanIndex', ...
        'ShotIndex', ...
        'PointIndex', ...
        'Time_s', ...
        'Voltage_mV', ...
        'dV_mV', ...
        'Counts_per_bin', ...
        'TemplateShape_counts_per_bin_for_A1', ...
        'FitCounts_per_bin', ...
        'Residual_counts_per_bin', ...
        'Afit_total_counts_above_baseline', ...
        'B0fit_counts_per_bin', ...
        'mu_ref_mV', ...
        'sigma_ref_mV', ...
        'f_ref', ...
        'tau_ref_mV'};

    C = headers;

    C(end+1,:) = { ...
        'metadata', ...
        char(dataSetType), ...
        char(exportTimestamp), ...
        char(string(transitionName)), ...
        char(string(libFile)), ...
        '', '', '', '', '', '', '', '', '', '', '', '', ...
        mu_ref, sigma_ref, f_ref, tau_ref};

    if ~isempty(exportRows)
        C = [C; exportRows];
    end

    writecell(C, fullpath);

    fprintf('Exported library-fit plot/fit data to:\n%s\n', fullpath);
end

function exportRows = append_library_singlefit_param_rows( ...
    exportRows, dataSetType, exportTimestamp, transitionName, libFile, ...
    scanIndex, shotIndex, Afit, mu_ref, sigma_ref, b0fit, f_ref, tau_ref)

    paramNames  = { ...
        'Afit_total_counts_above_baseline', ...
        'mu_ref_mV', ...
        'sigma_ref_mV', ...
        'b0fit_counts_per_bin', ...
        'f_ref', ...
        'tau_ref_mV'};

    paramValues = [Afit, mu_ref, sigma_ref, b0fit, f_ref, tau_ref];

    for k = 1:numel(paramNames)
        exportRows(end+1,:) = { ...
            'fit_parameter', ...
            char(dataSetType), ...
            char(exportTimestamp), ...
            char(string(transitionName)), ...
            char(string(libFile)), ...
            scanIndex, ...
            shotIndex, ...
            k, ...
            '', ...
            '', ...
            '', ...
            '', ...
            '', ...
            '', ...
            '', ...
            paramValues(k), ...
            '', ...
            '', ...
            '', ...
            '', ...
            paramNames{k}};
    end
end