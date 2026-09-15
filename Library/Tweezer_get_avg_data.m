function avgDataset = Tweezer_get_avg_data(analyVar,indivDataset)
    avgDataset = cell(1,length(analyVar.uniqScanList));
    for uniqScanIter = 1:length(analyVar.uniqScanList)
        % Create a cell containing the scanned variable for each batch (set of images)
        indVars = arrayfun(@(x) indivDataset{x}.imagevcoAtom,analyVar.posOccurUniqVar{uniqScanIter},'UniformOutput',0);
        % Save a vector of all the unique scan pnts for comparison
        avgDataset{uniqScanIter}.simScanIndVar = double(unique(int32(cell2mat(indVars)*analyVar.compPrec)))*1/analyVar.compPrec;
    end
end