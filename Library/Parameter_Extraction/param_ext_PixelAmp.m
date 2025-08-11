function [indivDataset] = param_ext_PixelAmp(analyVar,indivDataset)
% Function to extract the pixel amplitude from each 
% image found from the 2D fit parameters.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS:
%   indivDataset - Output will add the fields below to the structure
%                  indivDatset
%                     Pixel Amp - The amplitude in x and y axis for each
%                     image.

%%%%%%%-----------------------------------%%%%%%%%%%
% Cloud radius accounting for finite resolution of the imaging system
% Reference - Mi Yan's PhD thesis Appendix A (eq. A.6)
%effPixelSize = analyVar.sizefactor*1e6; %um, pixel size

%% Cloud parameters needed
%%%%%%%-----------------------------------%%%%%%%%%%

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
for basenameNum = 1:analyVar.numBasenamesAtom
    % Preallocate nested loop variables
    [cldAmpX cldAmpY]...
        = deal(NaN(sum(analyVar.LatticeAxesFit),indivDataset{basenameNum}.CounterAtom));
    
    % Processes all the image files in the current batch
    for k = 1:indivDataset{basenameNum}.CounterAtom;
%% Find Cloud Amp of each window
% Find size in each window (if multiple windows)
%%%%%%%-----------------------------------%%%%%%%%%%
        % Amp of Cloud
        cldAmpX(:,k) = indivDataset{basenameNum}.All_fitParams{k}{1}.Amp;
        cldAmpY(:,k) = indivDataset{basenameNum}.All_fitParams{k}{1}.Amp;
    end
    
    % Save cloud radius for each window into the indivDataset structure
    indivDataset{basenameNum}.cloudAmpX = cldAmpX;
    indivDataset{basenameNum}.cloudAmpY = cldAmpY;
end