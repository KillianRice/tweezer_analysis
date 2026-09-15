function InitCondCell = get_fit_params(analyVar,OD_Image,OD_Fit_ImageCell,indivDataset,basenameNum,k)
% Function to return the Initial conditions to be when fitting each individual
% cloud area. 
%
% INPUTS:
%   analyVar         - Structure containing all the global variables needed for analysis
%   OD_Image         - Processed OD image from the experiment (full image)
%   OD_Fit_ImageCell - Cell of OD images that will be fitted (may be
%                      smaller than matrix in OD_Image)
%   indivDataset     - Cell of structures containing all batch specific data
%   basenameNum      - Iteration variable to index data specific to the batch
%                      being analyzed
%   k                - Iteration variable to index data specific to the image
%                      being analyzed 
%
% OUTPUTS:
%   InitCondCell - Cell containing the initial guesses for non-linear
%                  fitting for each cloud image defined in roiWin_Index

% Updated 2026.8.21 To include a 1D integration fitting form for tweezer
% images
    
    %%  Initialize Variables
    InitCond = zeros(length(analyVar.InitCondList),1);
    
    %%  Smooth full OD image for guessing slopes and offset
    Zguess = conv2(OD_Image,fspecial('gaussian',analyVar.gaussFiltAmp,analyVar.gaussFiltSig),'same');
    
    %% ================================================================
    %  DIMENSIONAL-REDUCTION INITIAL GUESS
    %  ================================================================
    
    if isfield(analyVar,'dimenReduc1D') && ...
            analyVar.dimenReduc1D == 1
    
        InitCondCell = ...
            cell(size(OD_Fit_ImageCell));
    
        for n = 1:numel(OD_Fit_ImageCell)
    
            OD = ...
                OD_Fit_ImageCell{n};
    
            InitCondCell{n} = ...
                integration1DGuess( ...
                    analyVar, ...
                    OD, ...
                    InitCond, ...
                    indivDataset, ...
                    basenameNum, ...
                    k, ...
                    Zguess);
    
        end
    
        return;
    
    end
    
    %% ================================================================
    %  REGULAR 2D FITTING INITIAL GUESS (LEGACY)
    %  ================================================================
    
    %% Get Initial Conditions
    InitCondCell = cellfun(@(x) feval(str2func(analyVar.InitCase),analyVar,x,InitCond,indivDataset,basenameNum,k,Zguess)...
                            ,OD_Fit_ImageCell,'UniformOutput',0);
end


%% ================================================================
%  Initial Condition Functions Based on Form
%   Pure
%   Bimodal
%  ================================================================
function InitCond = Pure(analyVar,OD_Image,InitCond,indivDataset,basenameNum,k,Zguess)%#ok<DEFNU>
    %% Find initial guess for each window
    % Amp, xCntr, and yCntr
    % Find positions of variables in analyVar.findFromPeak in InitCondCell
    [~,cntrInd] = intersect(analyVar.InitCondList,analyVar.findFromPeak);
    % Assign peak and peak position to variables in findFromPeak
    [InitCond(cntrInd)] = funcPeakGuess(analyVar,OD_Image);
    
    % sigX and sigY (thermal portion widths)
    [InitCond([2 3])] = indivDataset{basenameNum}.sigParamAtom(k);
    
    % Linear background (Offset, SlopeX, and Slope Y)
    % Using entire ROI image not just cloud so these values should be the same for all windows
    [InitCond([6 7 8])] = linearBack(analyVar,Zguess);
end

function InitCond = Bimodal(analyVar,OD_Image,InitCond,indivDataset,basenameNum,k,Zguess)%#ok<DEFNU>
    %% Find inital guess for each window
    % Amp_BEC, xCntr, and yCntr
    % Find positions of variables in analyVar.findFromPeak in InitCondCell
    [~,cntrInd] = intersect(analyVar.InitCondList,analyVar.findFromPeak);
    % Assign peak and peak position to variables in findFromPeak
    [InitCond(cntrInd)] = funcPeakGuess(analyVar,OD_Image);
    
    % sigX_BEC and sigY_BEC
    [InitCond([2 3])] = size(OD_Image)/indivDataset{basenameNum}.sigBECParamAtom(k);
    
    % Amp (thermal amplitude)
    InitCond(4)       = InitCond(1)*analyVar.ampBimodalGuess;
    
    % sigX and sigY (thermal portion widths)
    [InitCond([6 5])] = indivDataset{basenameNum}.sigParamAtom(k);
    
    % Linear background (Offset, SlopeX, and Slope Y)
    % Using entire ROI image not just cloud so these values should be the same for all windows
    [InitCond([9 10 11])] = linearBack(analyVar,Zguess);
end



%% ================================================================
%  Common functions between Bimodal and Single routines
%  ================================================================
function peakGuess = funcPeakGuess(analyVar,OD_Image)
    %% Find peak amplitude and position 
    %Cloud center guess by filtering with guassian (same filter as in Zguess)
    try
        cntr = FastPeakFind(OD_Image',0,fspecial('gaussian',analyVar.gaussFiltAmp,analyVar.gaussFiltSig),2);
    catch %#ok<CTCH>
        % Sometimes FastPeakFind fails to find a peak (when analyzing windows
        % with no population) so default to middle of the window.
        cntr = round(size(OD_Image)'./2);
    end
    % Assign highest amplitude to the first variable and x/y centers accordingly
    peakGuess = [OD_Image(cntr(1),cntr(2)),cntr(1),cntr(2)];
end

function linearGuess = linearBack(analyVar,Zguess)
    %% Average data in corners to find overal linear background
    
    % Find average value in the corners (size determined by analyVar.NoiseNumVec)
    averageYLowXLow   = sum(sum(Zguess(1:analyVar.NoiseNumVec,1:analyVar.NoiseNumVec)))/analyVar.NoiseNumVec^2;
    averageYLowXHigh  = sum(sum(Zguess(1:analyVar.NoiseNumVec,...
        size(Zguess,2) - analyVar.NoiseNumVec:size(Zguess,2))))/analyVar.NoiseNumVec^2;
    averageYHighXLow  = sum(sum(Zguess(size(Zguess,1) - analyVar.NoiseNumVec:size(Zguess,1),...
        1:analyVar.NoiseNumVec)))/analyVar.NoiseNumVec^2;
    averageYHighXHigh = sum(sum(Zguess(size(Zguess,1) - analyVar.NoiseNumVec:size(Zguess,1),...
        size(Zguess,2) - analyVar.NoiseNumVec:size(Zguess,2))))/analyVar.NoiseNumVec^2;
    
    % Estimate image slope from averaged element values in the corners of the image
    % Used in fitting to compensate for image not being perpendicular to camera.
    linearGuess = [(averageYLowXLow + averageYLowXHigh + averageYHighXLow + averageYHighXHigh)/4
                  ((averageYLowXHigh + averageYHighXHigh) - (averageYLowXLow + averageYHighXLow))/(2*size(Zguess,2))
                  ((averageYHighXLow + averageYHighXHigh) - (averageYLowXLow + averageYLowXHigh))/(2*size(Zguess,1))];
end


function InitCond = integration1DGuess( ...
    analyVar,OD_Image,InitCond, ...
    indivDataset,basenameNum,k,Zguess)
    
    %% ---------------------------------------------------------------
    %  Integrated profiles
    
    xProfile = sum(OD_Image,1);
    yProfile = sum(OD_Image,2);
    
    %% ---------------------------------------------------------------
    %  Estimate cloud center
    
    [~,xCntr] = max(xProfile);
    [~,yCntr] = max(yProfile);
    
    %% ---------------------------------------------------------------
    %  Estimate amplitude
    % The integrated profile amplitude is larger than the original
    % 2D peak amplitude, so use the 2D peak as the amplitude guess.
    
    peakValue = ...
        max(OD_Image(:));
    
    %% ---------------------------------------------------------------
    %  Gaussian widths
    % Preserve the existing experimentally determined widths when
    % available. This is the safest initial guess and avoids changing
    % the behavior of previous analysis modes.

    sigX = 5;% size(OD_Image,2)/2;
    sigY = 5;% size(OD_Image,1)/2;
    
    %% ---------------------------------------------------------------
    %  Assign according to InitCondList

    parameterNames = ...
        analyVar.InitCondList;
    
    for p = 1:numel(parameterNames)
    
        switch lower(parameterNames{p})
    
            case {'amp','amp_bec'}
                InitCond(p) = peakValue;
    
            case {'sigx','sigx_bec'}
                InitCond(p) = sigX;
    
            case {'sigy','sigy_bec'}
                InitCond(p) = sigY;
    
            case {'xcntr','xcenter','x0'}
                InitCond(p) = xCntr;
    
            case {'ycntr','ycenter','y0'}
                InitCond(p) = yCntr;
    
        end
    
    end
    
    %% ---------------------------------------------------------------
    %  Background
    
    if strcmpi(analyVar.InitCase,'Pure')
    
        [InitCond([6 7 8])] = [0,0,0];%...
            %linearBack(analyVar,Zguess);
    
    elseif strcmpi(analyVar.InitCase,'Bimodal')
    
        [InitCond([9 10 11])] = ...
            linearBack(analyVar,Zguess);
    
    end
    
end