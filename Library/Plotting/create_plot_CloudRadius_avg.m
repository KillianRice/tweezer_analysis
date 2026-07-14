function create_plot_CloudRadius_avg(analyVar,avgDataset)
% Function to plot the cloud radius
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS:
%   Creates plots showing the cloud radius

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar    = analyVar.funcDataScale(avgDataset.imagevcoAtom);
    CloudRadX = avgDataset.cloudRadX(1,:);
    CloudRadY = avgDataset.cloudRadY(1,:);
    
%% Plot of radius in X & Y
%%%%%%%-----------------------------------%%%%%%%%%%
    figNum = analyVar.figNum.atomSize; 
    radLabel = {'Cloud Radius [um]', 'Cloud Radius [um]'}; 
    radTitle = {'', 'X - axis', 'Y - axis'};
    default_plot(analyVar,[1 analyVar.numBasenamesAtom],...
        figNum,radLabel,radTitle,analyVar.timevectorAtom,...
        repmat(indVar,1,2)',[CloudRadX; CloudRadY]);
    if analyVar.fitModel == 'PureGaussian'
        str={'Fit Model: $$Z= A(\exp^{-(X^2/2\sigma_x^2-Y^2/2\sigma_y^2)})$$'};
        annotation('textbox','interpreter','latex','String',str,'FitBoxToText','on')
    end
end