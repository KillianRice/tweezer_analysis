% Script to run complete analysis of Neutral table data
% 
% Please be sure to have the correct date defined in AnalysisVariables
% prior to running.

% clear current figures
close all

% Load variables and file data
analyVar = AnalysisVariables;
if analyVar.numBasenamesAtom > 5
    disp(['About to analyze ' num2str(analyVar.numBasenamesAtom) ' scans, are you sure you want to continue?']);
    s = input('Y/N >','s');
    while s ~= 'Y' && s ~= 'N' && s ~= 'y' && s ~= 'n'
        s = input('Y/N >','s');
        disp(s)
    end
    if s == 'n' || s == 'N'
        return;
    end
end
indivDataset = get_indiv_batch_data(analyVar);

if analyVar.UseTweezer == 1
    % 1. make/save individual OD images and Raw ROI cuts
    imagefit_Backgrounds_PCA_V2(analyVar, indivDataset);

    % 2. Grab OD Total Counts for Histograms
    indivDataset = param_ext_ODTotalCounts(analyVar,indivDataset);
    
    % 3. Average matching scan points across scans
    avgDataset = imagefit_BuildAveragedScans(analyVar, indivDataset);
    
    % 4. Fit averaged OD images only (not really working for small images)
    if analyVar.fitODImage
        avgDataset = imagefit_NumDistFit_Averaged(analyVar, indivDataset, avgDataset);
    end
    
    % 5. Evaluate averaged fits
    if analyVar.SavePlotData == 1
        PlotData = imagefit_ParamEval_Averaged(analyVar, indivDataset, avgDataset);
    else
        imagefit_ParamEval_Averaged(analyVar, indivDataset, avgDataset);
    end
else

    % Background fitting
    imagefit_Backgrounds_PCA(analyVar,indivDataset)
    
    % Functional fitting
    imagefit_NumDistFit(analyVar,indivDataset)
    
    % Plotting routine
    if analyVar.SavePlotData == 1
        % Save data from output
        PlotData = imagefit_ParamEval(analyVar,indivDataset);
    else
        imagefit_ParamEval(analyVar,indivDataset);
    end 
end