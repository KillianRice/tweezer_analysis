function create_plot_COM(analyVar,indivDataset)
% Function to plot the cloud center of mass
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS: plots center of mass 
% Created by SKK Mar 31, 2025.

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
for basenameNum = 1:analyVar.numBasenamesAtom
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar  = analyVar.funcDataScale(indivDataset{basenameNum}.imagevcoAtom);
    X_com   = indivDataset{basenameNum}.x_com(1,:);
    Y_com   = indivDataset{basenameNum}.y_com(1,:);
   
%% Plot of radius in X & Y
%%%%%%%-----------------------------------%%%%%%%%%%
    figNum = analyVar.figNum.COM; 
    radLabel = {'COM [m]', 'COM [m]'};
    radTitle = {'', 'X - axis', 'Y - axis'};
    default_plot(analyVar,[basenameNum analyVar.numBasenamesAtom],...
        figNum,radLabel,radTitle,analyVar.timevectorAtom,...
        repmat(indVar,1,2)',[X_com; Y_com]);
     str={'Effective Pixelsize',analyVar.sizefactor};
     annotation('textbox','interpreter','latex','String',str,'FitBoxToText','on')
end
end