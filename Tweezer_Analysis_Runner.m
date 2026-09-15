% Script to run Tweezer image analyses.
% 
%% clear current figures
close all

%% Load variables and file data
analyVar = TweezerAnalysisVariables;
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
indivDataset = Tweezer_get_indiv_batch_data(analyVar);
avgDataset = Tweezer_get_avg_data(analyVar,indivDataset);
Tweezer_imagefit_ParamEval(analyVar,indivDataset,avgDataset);