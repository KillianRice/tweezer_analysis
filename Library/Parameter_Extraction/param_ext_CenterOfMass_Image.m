function [indivDataset] = param_ext_CenterOfMass_Image(analyVar,indivDataset)
%This program is designed to read OD files saved after
%imagefit_Analysis_Runner is successfully run, which saves the OD files for
%each batch file.
%It will return center of mass of [x_cm, y_cm].
%
% INPUTS: AnalyVar, indivDataset
%
% OUTPUTS:
%   [x_cm, y_cm]
% Created by SKK Mar 31, 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through each batch file listed in basenamevectorAtom
for basenameNum = 1:analyVar.numBasenamesAtom
    % Preallocate nested loop variables
    [x_com, y_com]...
        = deal(NaN(1,indivDataset{basenameNum}.CounterAtom));
    
    % Process all the data files in this batch
    for k = 1:indivDataset{basenameNum}.CounterAtom;
%% Retrieve OD image from file
%%%%%%%-----------------------------------%%%%%%%%%%
        if exist([analyVar.analyOutDir char(indivDataset{basenameNum}.fileAtom(k)) analyVar.ODimageFilename],'file')
            OD_Image_Single = dlmread([analyVar.analyOutDir char(indivDataset{basenameNum}.fileAtom(k)) analyVar.ODimageFilename]);
        else
            % Error if background subtraction not run (means no OD file saved)
            callStruct = dbstack(0);
            error('%s Cannot find \n\t%s',...
                callStruct.name,[indivDataset{basenameNum}.fileAtom{k} analyVar.ODimageFilename]);
        end
        [x_com(:,k), y_com(:,k)] = ReturnCOM_ODImage(OD_Image_Single);
    end % end loop through each dataset
    indivDataset{basenameNum}.x_com = x_com.*analyVar.pixelsize;% SI units real units
    indivDataset{basenameNum}.y_com = y_com.*analyVar.pixelsize;% SI units real units
end     % end loop through master batch file
