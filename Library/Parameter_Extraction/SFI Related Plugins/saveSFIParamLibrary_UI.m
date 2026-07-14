function saveSFIParamLibrary_UI(indivDataset, analyVar)
% Save the GLOBAL summed SFI fit parameters to the shared TSV library file.
%
% Legacy schema:
%   timestamp, scan_i, shot_j, transition_name, A, mu, sigma, b0, f, tau, dt
%
% For the new global workflow:
% - scan_i is written as NaN
% - shot_j is written as NaN
% - dt stores globalSummedSFISpectra_dV if available

    %#ok<INUSD> analyVar kept for compatibility with existing calls

    % ------------------------------------------------------------
    % Find the first dataset entry containing the global summed fit
    % ------------------------------------------------------------
    idx = [];
    for i = 1:numel(indivDataset)
        if isempty(indivDataset{i}) || ~isstruct(indivDataset{i})
            continue;
        end

        if isfield(indivDataset{i}, 'globalSkewFitParams_summed') && ...
                ~isempty(indivDataset{i}.globalSkewFitParams_summed)
            idx = i;
            break;
        end
    end

    if isempty(idx)
        tmpFig = uifigure('Visible','off');
        uialert(tmpFig, ...
            'No global summed fit parameters were found to save.', ...
            'Nothing to Save');
        return;
    end

    % ------------------------------------------------------------
    % Ask user for parameter-set / transition name
    % ------------------------------------------------------------
    answer = inputdlg( ...
        {'Enter parameter set / transition name:'}, ...
        'Save SFI Fit Parameters', ...
        [1 60]);

    if isempty(answer)
        return;
    end

    transitionName = strtrim(answer{1});
    if isempty(transitionName)
        tmpFig = uifigure('Visible','off');
        uialert(tmpFig, ...
            'Parameter set name cannot be empty.', ...
            'Missing Name');
        return;
    end

    % ------------------------------------------------------------
    % Pull fit parameters and metadata
    % ------------------------------------------------------------
    p = indivDataset{idx}.globalSkewFitParams_summed;

    if numel(p) < 6
        error('globalSkewFitParams_summed is malformed. Expected at least 6 parameters.');
    end

    A   = p(1);
    mu  = p(2);
    sig = p(3);
    b0  = p(4);
    f   = p(5);
    tau = p(6);

    % Legacy schema placeholders for global fits
    scan_i = NaN;
    shot_j = NaN;

    if isfield(indivDataset{idx}, 'globalSummedSFISpectra_dV') && ...
            ~isempty(indivDataset{idx}.globalSummedSFISpectra_dV)
        dt = indivDataset{idx}.globalSummedSFISpectra_dV;
    else
        dt = NaN;
    end

    if isfield(indivDataset{idx}, 'globalSummedSFISpectraContributors') && ...
            ~isempty(indivDataset{idx}.globalSummedSFISpectraContributors)
        nContrib = size(indivDataset{idx}.globalSummedSFISpectraContributors, 1);
    else
        nContrib = NaN;
    end

    % ------------------------------------------------------------
    % Write / append to library file
    % ------------------------------------------------------------
    here = fileparts(mfilename('fullpath'));
    outFile = fullfile(here, 'sfi_fit_parameter_library.txt');

    legacyHeader = sprintf(['timestamp\tscan_i\tshot_j\ttransition_name\tA\tmu\t' ...
                            'sigma\tb0\tf\ttau\tdt\n']);

    expectedHeader = ['timestamp' char(9) ...
                      'scan_i' char(9) ...
                      'shot_j' char(9) ...
                      'transition_name' char(9) ...
                      'A' char(9) ...
                      'mu' char(9) ...
                      'sigma' char(9) ...
                      'b0' char(9) ...
                      'f' char(9) ...
                      'tau' char(9) ...
                      'dt'];

    % If file does not exist, create it with legacy header
    if ~exist(outFile, 'file')
        fid = fopen(outFile, 'w');
        if fid < 0
            error('Could not open output file for writing: %s', outFile);
        end
        fprintf(fid, legacyHeader);
        fclose(fid);
    else
        % Sanity check: warn if header is unexpected
        fid = fopen(outFile, 'r');
        if fid >= 0
            firstLine = fgetl(fid);
            fclose(fid);

            if ischar(firstLine) && ~strcmp(strtrim(firstLine), expectedHeader)
                warning('Library file header does not match expected legacy format. Appending anyway.');
            end
        end
    end

    fid = fopen(outFile, 'a');
    if fid < 0
        error('Could not open output file for appending: %s', outFile);
    end

    ts = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    fprintf(fid, '%s\t%.0f\t%.0f\t%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n', ...
        ts, scan_i, shot_j, transitionName, A, mu, sig, b0, f, tau, dt);

    fclose(fid);

    % ------------------------------------------------------------
    % Confirmation
    % ------------------------------------------------------------
    msg = sprintf(['Saved parameter set:\n\n' ...
                   '  Name: %s\n' ...
                   '  File: %s\n' ...
                   '  Contributors: %g'], ...
                   transitionName, outFile, nContrib);

    tmpFig = uifigure;
    uialert(tmpFig, msg, 'Save Complete');
end