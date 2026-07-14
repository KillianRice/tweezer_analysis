function evolAxH = create_plot_evol_avg(analyVar,avgDataset)
% Function to plot the OD images of a scan on a subplot.
% This function does no parameter extraction from the fit but serves only
% to show the evolution of the cloud through a scan.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS:
%   evolAxH - Cell of vectors containing handles to the subplot axes. 
%             Used for setting standard color limits on pcolor plots.
%

%% Preallocate loop variables
evolAxH = zeros(1,avgDataset.CounterAtom);

%% Loop through each avg batch
for k = 1:avgDataset.CounterAtom
    %% Preallocate loop variables
    
        % Processes all the image files in this batch/scan
       figure(analyVar.figNum.atomEvol);
            %% Plot Evolution
            evolAxH(k) = subplot(avgDataset.SubPlotRows,avgDataset.SubPlotCols,k);
            pcolor(avgDataset.All_OD_Image{k});
            colorbar;
			shading flat;
			axis equal tight;
            %%% Plot axis details
            title(strcat(num2str(avgDataset.imagevcoAtom(ceil(k/avgDataset.numTweezers))),[' ' analyVar.xDataUnit]));
    
end
hold on; grid on; axis on
set(gcf,'Name',['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(1))]);
mtit(['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(1))],'FontSize',16,'zoff',.05,'xoff',-.01)