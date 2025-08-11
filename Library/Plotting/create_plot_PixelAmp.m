function create_plot_PixelAmp(analyVar,indivDataset)
% Function to plot the cloud radius
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS:
%   Creates plots showing the cloud amplitude in pixels

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
for basenameNum = 1:analyVar.numBasenamesAtom
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar    = analyVar.funcDataScale(indivDataset{basenameNum}.imagevcoAtom);
    CloudAmpX = indivDataset{basenameNum}.cloudAmpX(1,:);
    CloudAmpY = indivDataset{basenameNum}.cloudAmpY(1,:);
    
%% Plot of radius in X & Y
%%%%%%%-----------------------------------%%%%%%%%%%
    figNum = analyVar.figNum.Amp; 
    ampLabel = {'Cloud Amp', 'Cloud Amp'}; 
    ampTitle = {'', 'X - axis', 'Y - axis'};
    default_plot(analyVar,[basenameNum analyVar.numBasenamesAtom],...
        figNum,ampLabel,ampTitle,analyVar.timevectorAtom,...
        repmat(indVar,1,2)',[CloudAmpX; CloudAmpY]);
    if analyVar.fitModel == 'PureGaussian'
        str={'Fit Model: $$Z= A(\exp^{-(X^2/2\sigma_x^2-Y^2/2\sigma_y^2)})$$'};
        annotation('textbox','interpreter','latex','String',str,'FitBoxToText','on')
    end
end
end