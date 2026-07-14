% =====================================================================
% Helper: read library file and pick most recent row for transitionName
% =====================================================================
function [mu_ref, sigma_ref, alpha_ref] = local_get_template(libFile, transitionName)
    if ~exist(libFile,'file')
        error('Library file not found: %s', libFile);
    end

    opts = detectImportOptions(libFile,'FileType','text');
    opts.Delimiter = '\t';
    T = readtable(libFile, opts);

    required = ["transition_name","mu","sigma","alpha"];
    if ~all(ismember(required, string(T.Properties.VariableNames)))
        error('Library file missing required columns: %s', strjoin(required, ', '));
    end

    names = string(T.transition_name);
    target = string(transitionName);

    rows = find(names == target);
    if isempty(rows)
        error('Transition name not found in library: %s', target);
    end

    r = rows(end); % most recent match

    mu_ref    = T.mu(r);
    sigma_ref = T.sigma(r);
    alpha_ref = T.alpha(r);

    if ~isfinite(mu_ref) || ~isfinite(sigma_ref) || ~isfinite(alpha_ref) || sigma_ref <= 0
        error('Invalid template values in library for "%s": mu=%g sigma=%g alpha=%g', target, mu_ref, sigma_ref, alpha_ref);
    end
end