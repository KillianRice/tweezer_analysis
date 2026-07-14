function tf = local_has_fitted_sfi_integral(indivDataset)
% Returns true if it looks like the new SFI fitting pipeline has already populated
% per-shot fitted integrals (Afit) and related fields.

    tf = false;
    if isempty(indivDataset), return; end

    for i = 1:numel(indivDataset)
        indiv = indivDataset{i};
        if isempty(indiv) || ~isstruct(indiv)
            continue;
        end

        % Strong indicator: library/raw all-shots fit was run
        if isfield(indiv,'sfiFitA') && ~isempty(indiv.sfiFitA) && any(isfinite(indiv.sfiFitA))
            tf = true;
            return;
        end

        % Alternate indicator: fit params stored per shot
        if isfield(indiv,'sfiFitParams') && ~isempty(indiv.sfiFitParams)
            % If any cell entry is non-empty, treat as already fitted
            try
                if any(~cellfun(@isempty, indiv.sfiFitParams))
                    tf = true;
                    return;
                end
            catch
            end
        end

        % Weak indicator: sfiIntegral already computed (cell or numeric) and non-empty
        if isfield(indiv,'sfiIntegral') && ~isempty(indiv.sfiIntegral)
            if isnumeric(indiv.sfiIntegral) && any(isfinite(indiv.sfiIntegral(:)))
                tf = true; return;
            end
            if iscell(indiv.sfiIntegral)
                try
                    vals = cellfun(@(x) double(x), indiv.sfiIntegral, 'UniformOutput', false);
                    vals = [vals{:}];
                    if any(isfinite(vals))
                        tf = true; return;
                    end
                catch
                    % If it's a cell but not numeric, ignore
                end
            end
        end
    end
end