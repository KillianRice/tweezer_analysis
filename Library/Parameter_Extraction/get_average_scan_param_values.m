function paramVals = get_average_scan_param_values(analyVar, indivDataset, basenameNum)
% Function will grab the parameter that was scanned over (such as
% imagevcoatom) or the ID of the scan (if doing dummy scans) so that
% imagefit_BuildAverageScans can collect the data in the intended manner


%% Order by scan parameter -- imagevcoAtom
fieldName = 'imagevcoAtom';% analyVar.avgScanParamField;

if isfield(indivDataset{basenameNum}, fieldName)
    paramVals = indivDataset{basenameNum}.(fieldName);
else
    error('Could not find scan parameter field "%s" in indivDataset{%d}.', fieldName, basenameNum);
end

paramVals = paramVals(:);

%% Order by ID of file
%%images are averaged over the same file rather than different files
if analyVar.dummyScan
    for idx = 1:length(paramVals)
        paramVals(idx) = analyVar.meanListVar(basenameNum);  %%The id of the scan
    end
end


%% Check length of params to number of scans in a file
if numel(paramVals) ~= indivDataset{basenameNum}.CounterAtom
    error('Parameter vector length does not match CounterAtom for basenameNum %d.', basenameNum);
end

end