function analyVar = TweezerAnalysisVariables()
% Master file to control all similar variables to be passed between
% background analysis, cloud fitting, and graphing routines.
%
% INPUTS:
%   none
%
% OUTPUTS:
%   analyVar - Structure containing all the variables defined in
%              this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS TO EDIT for each experiment.
CameraType = 1;                 % set to 1 to use Zyla4.2 sideview camera and 0 to use the PixelFly
quantumNumberN = 44;            % principal quantum number n when doing Rydberg experiments.
state = '1S0';                  % term symbol for rydberg state (2S+1)L(J)
NeutExpDir      = 'Raw_Data';
analyPrefix     = '_twzrAccStudyANDAIDecay';                                % This is the suffix name of analysis folder.
analyOutputName = 'Analysis';
UseMCS = 0; 
UseWavemeter = 0;
UseImages = 1;
SavePlotData = 1;               % Keep this 1

%% Global Exp. Settings
isotope = 88;                                                               % This is fixed for Tweezer experiment.
detuning = 0;                                                               % s^-1, image beam detuning (as of 7/1/15).

%% Rydberg properties
quantumDefect = 0;
switch state
    case '3S1'
        quantumDefect = 3.371;
    case '1S0'
        quantumDefect = 3.26896;
    case '1D2'
        quantumDefect = 2.3807;
    case '3D1'
        quantumDefect = 2.658;
    case '3D2'
        quantumDefect = 2.636;
    case '3D3'
        quantumDefect = 2.63;
end
nStar = quantumNumberN - quantumDefect;

%% Camera Settings
if CameraType == 0
    binHorizontal  = 2;%binning done by camera when taking images
    binVertical    = 2;
    matrixSize     = [1280/binVertical 1024/binHorizontal]; % Matrix size of camera output: Set this to be the same as PixelFly dimensions.
    CameraRes  = 15; %um
    pixelOnCam = 6.7*10^(-6); %m
    MagImgSystem = 1;
    bin = binHorizontal;
    pixelsize  = bin*pixelOnCam/MagImgSystem; %m/px
end
if CameraType == 1
    binHorizontal  = 2;     % binning done by camera when taking images
    binVertical    = 2;
    matrixSize     = [400/binVertical 400/binHorizontal];   % Matrix size of camera output: Set this to be the same as Zyla dimensions.
    CameraRes  = 0.8;       % um measured with free falls measurements.
    pixelOnCam = 6.5*10^(-6);   % m
    MagImgSystem = 8;       % 25 mm in-vacuo & 200 mm for tubelens. Change Mag when using different lens, like 200 mm for Tweezers.
    bin = binHorizontal;
    pixelsize  = bin*pixelOnCam/MagImgSystem;   % m/px
end 

%% LINESHAPE FITTING
% Called within imagefit_ParamEval
% Allows secondary analysis of parameters extracted from the number distribution fit
% Functions ahould expect input arguments of analyvar, indivDataset, and avgDataset from imagefit_ParamEval
%
% Functions are selected by the booleans in lcl_logicFitLine (i.e. [1 0 0 0] triggers lcl_validFitLine(1))
%   Spectrum_Fit - Fits spectrum to gaussian (or lorentzian) and estimates rabi freq. based on width
%   KapitzaDirac - Fits number oscillations in 2hk peaks to calibrate lattice depth
%   Cloud_Pos    - Plot spatial center of cloud (option to fit oscillation for trap freq. measurements)
% Future Plans
%   lifetime_n   - (n = 1,2,3) Fits number decay to decaying exponential with n-body decay constant
lcl_validFitLine = {'CheckRawImages_v2',...                  % 01 
                    'CheckTweezersROI',...
                    'TweezerROILineIntegrals',...
                    'FitTweezerLineIntegrals',...
                    };

plugInVec = [2,3];

%% Check if roifile exists in /out
roiFile = 'out/ROI_Settings_current.mat';
if exist(roiFile, 'file')
    % Load ROI settings
    roiData = load(roiFile);
    roiCenters = roiData.roiSettings.centers;
    roiRadius  = roiData.roiSettings.radius;
else
    % Ask user whether to skip ROI loading
    answer = questdlg( ...
        'ROI settings file was not found. Skip ROI loading?', ...
        'ROI Settings Not Found', ...
        'Yes', 'No', 'Yes');
    if strcmp(answer, 'Yes')
        fprintf('Skipping ROI loading.\n');
        roiCenters = [];
        roiRadius = [];
    else
        error('ROI settings file not found: %s', roiFile);
    end
end

%% Common Plotting flags
    lcl_logicFitLine = zeros(1,length(lcl_validFitLine)); 
if isempty(plugInVec )~= 1
    lcl_logicFitLine(plugInVec) = 1;
end

%% Check for valid lineshape fit arguments
fitLineFunc = lcl_validFitLine(nonzeros(lcl_logicFitLine.*(1:length(lcl_validFitLine))));
if  isempty(intersect(lcl_validFitLine,fitLineFunc)) && ~isempty(fitLineFunc)
    error('Invalid option in Fitlineshape.')
end

%% PicoScope settings (if using)
plotCounts = 0; % look for scope traces from picoscope
SumCounts = 0;

%% SR400 Photon Counter (if using)
plotCounts_SR400 = 0; %photon counter

%% What kind of files are expected?
%%-----------------------------------------------------------------------%%
dataAtom = char('atoms.bny');
dataBack = char('back.bny');
dataMCS = char('_mcs_counts.mcs');
%%-----------------------------------------------------------------------%%

%% Define Location of data and analysis
% Define default folder names for directory heirarchy
lcl_analyDir = pwd; % save Analysis Folder location (USE THIS NORMALLY)
disp(pwd);

% Determine directory where all raw data files are saved (expected to mirror folder structure of Analysis folder)
dataDir  = [strrep(strrep(lcl_analyDir,analyPrefix,''),[filesep 'Analysis' filesep],[filesep 'Raw_Data' filesep]) filesep];
if not(exist(dataDir,'dir'))
    error('\nData Directory: %s\n not found. Please check analysis directory and/or analyPrefix and try again.',dataDir)
end
dataDirName = regexp(dataDir,filesep,'split');
dataDirName = regexp(dataDirName{end - 1},'_','split');
dataDirName = dataDirName{1};
% Add the library folders to the path
addpath(genpath([pwd filesep 'Library']));
rmpath([pwd filesep 'Library' filesep 'Archive']);

%% LABVIEW BATCHFILE VARIABLES
%%-----------------------------------------------------------------------%%
% Variables defined in lines of each dataset in the master batch file.
lcl_masterBatchAtomVar = {
    'basenamevectorAtom'    % 01 name of data set
    'timevectorAtom'        % 02 four digit time stamp
    'ScanIDVarAtom'         % 03 ID used for averaging
    'redPower'              % 04
    'uvPower'               % 05
    'expTime'               % 06
    'droptimeAtom'          % 07 s, drop time before imaging
    'roiWinRadAtom'         % 08 pixels, radius of atoms
    'cloudWinRadAtom'       % 09 pixels, radious of image to crop (with background)
    'cloudColCntrAtom'      % 10 pixel, horizontal center of image to crop
    'cloudRowCntrAtom'      % 11 pixel, veritcal center of image to crop
    'mcs_roiStart'          % 12 first bin of MCS roi
    'mcs_roiEnd'            % 13 last bin of MCS roi
    };

colHeadersAtom = lcl_masterBatchAtomVar;
lcl_masterBatchAtomVar = lcl_masterBatchAtomVar';
lcl_masterBatchBackVar = {'basenamevectorBack' 'unusedBack'};

% Variables defined in lines of the individual batch files
indivBatchAtomVar = {
    'fileAtom' %                    01  name of the image file
    'imagevcoAtom' %                02  corresponding value of the independent parameter
    'principalQuantumAtom' %        03  principle quantum number n
    'angularQuantumAtom' %          04  principle quantum number \el
    'ODTHold' %                     05  odt hold time
    'RampDelay' %                   06  rydberg hold time for lifetime measurements
    'VCA_1_Voltage' %               07  3 beam ODT vca 1 static voltage; control relatice power between 3 beams
    'VCA_2_Voltage' %               08  3 beam ODT VCA 2 static voltage; control absolute power of all 3 beams
    'initialTrapDepthAtom' %        09  initial trap depth in volts before evaporation
    'TrapPower' %                   10  final trap depth in volts after evaporation
    'synthFreq' %                   11  synth frequency driving 640nm cat's eye aom
    'numberAtom' %                  12  number of atoms measured by labview
    'tempXAtom' %                   13  TOF x temperature measured by labview
    'tempYAtom' %                   14  TOF y temperature measured by labview
    'fugacityguessAtom' %           15  left constant
    'sigParamAtom' %                16  left constant
    'sigBECParamAtom' %             17  left constant
    'WeightedBECPeakAtom' %         18  left constant
    'BECamplitudeParameterAtom' %   19  left constant
    'wavemeterAtom'             %   20  wavemeter reading
    'timestampAtom'            %    21  timestamp of image (YYYY.MM.DD - HH:MM:SS)
    };
indivBatchBackVar = {
    'fileBack'
    'imagevcoBack'
    'principalQuantumBack' %principle quantum number n
    'angularQuantumBack' %principle quantum number \el
    'ODTHoldTimeBack' %odt hold time
    'RydbergHoldTimeBack' %rydberg hold time for lifetime measurements
    'VCA_1_VoltageBack' %3 beam odt vca 1 static voltage
    'VCA_2_VoltageBack' %3 beam ODT VCA 2 static voltage
    'initialTrapDepthBack' %initial trap depth in volts before evaporation
    'finalTrapDepthBack' %final trap depth in volts after evaporation
    'uvSynthFreqBack' %synth frequency driving 640nm cat's eye aom
    'numberBack' %number of atoms measured by labview
    'tempXBack' %TOF x temperature measured by labview
    'tempYBack' %TOF y temperature measured by labview
    'fugacityguessBack' %left constant
    'sigParamBack' %left constant
    'sigBECParamBack' %left constant
    'WeightedBECPeakBack' %left constant
    'BECamplitudeParameterBack' %left constant
    'wavemeterBack'             % wavemeter reading
    'timestampBack'            %    21  timestamp of image (YYYY.MM.DD - HH:MM:SS)
    };
% Grabs all the file names
basenamelistAtom = [dataDir 'Files_' strrep(dataDirName,'.','') '.txt'];
basenamelistBack = [dataDir 'Files_' strrep(dataDirName,'.','') '_Bg' '.txt'];

% Read in list of variable names and determine the format string for textscan
lcl_varFormatStr = cell(1,length(lcl_masterBatchAtomVar) );     % String is as long as number of variables
lcl_varFormatStr(1) = {'%s'};                                   % Define fixed variables
lcl_varFormatStr(cellfun('isempty',lcl_varFormatStr)) = {'%f'}; % Assume all other variables are floating point numbers
lcl_varFormatStr = horzcat(lcl_varFormatStr{:});                % Concatenate cells into single string for textscan

% Read in variables for each data set defined in the master batch file
lcl_masterBatchAtomData = textscan(fopen(basenamelistAtom),lcl_varFormatStr,'commentstyle','%'); %data set batch files for atoms
lcl_masterBatchBackData = textscan(fopen(basenamelistBack),'%s%f','commentstyle','%');         %data set batch files for backgrounds

%% Setup fields for averaging datasets together
%%-----------------------------------------------------------------------%%
% Used in imagefit_Plotting routine
% Allows averaging across multiple datasets with same parameters

meanListVar  = lcl_masterBatchAtomData{strcmpi(lcl_masterBatchAtomVar,'ScanIDVarAtom')}; % Variable used to identify similar scans
uniqScanList = unique(meanListVar,'stable'); % Unique values between all scans (maintains order of appearance in meanListVar)
posOccurUniqVar = arrayfun(@(x) find(meanListVar == x),uniqScanList,'UniformOutput',0);

%% PHYSICAL QUANTITIES
%%-----------------------------------------------------------------------%%
mass            = isotope*1.672621777e-27;                                  % kg, strontium mass
kBoltz          = 1.3806488e-23;                                            % J K^-1, Boltzmann's Constant
hbar            = 1.054571726e-34;                                          % J s, Reduced Planck's Constant
lambda          = 461e-9;                                                   % m, Sr 1S0->1P1 wavelength
lambdaLat       = 532e-9;                                                   % m, Tweezer wavelength. 
BohrRadius      = 5.2917721092e-11;                                         % m, Bohr Radius in meter
NaturalWidth    = 2*pi*30.5e6;                                              % s^-1, 1S0->1P1 FWHM
CrossSection    = 3*lambda^2/(2*pi);                                        % m^2, atom-photon cross-section, page 138 Foot Book.
AbsCross        = CrossSection*1/(1 + (2*detuning/NaturalWidth)^2);         % Pascal's PhD thesis (Eq. B.2)
RecoilEnergy    = hbar^2/(2*mass)*(2*pi/lambdaLat)^2;                       % One photon recoil energy of the lattice wavelength
Gravity         = 9.8;                                                      % m/s^2, acceleration due to gravity
Epsilon0        = 8.854187817e-12;                                          % F/m, permittivity of free space
SpeedOfLight    = 2.99792458e8;                                             % m/s, speed of light
aBohr           = 5.2917721067e-11;                                         % m, Bohr radius
QDefect         = 3.372;                                                    % unitless, Quantum Defect for 3S1 states, n>20
%%-----------------------------------------------------------------------%%
compPrec = 1e6; % Will round numbers to the 6th decimal place

%% This has to be at the end of the script.
%%% Enumerate number of Basenames
numBasenamesAtom = size(lcl_masterBatchAtomData{1},1);
% First need to match the batch filename variables to their values
lcl_masterBatchVars = cell2struct(cat(2,lcl_masterBatchAtomData,lcl_masterBatchBackData),cat(2,lcl_masterBatchAtomVar,lcl_masterBatchBackVar),2);
allVar   = who(); % analyVar is created with variables not starting with lcl_ (these are local variables that are not needed later in the analysis)
analyVar = catstruct(v2struct(cat(1,'fieldNames',allVar(cellfun('isempty',regexp(allVar,'\<lcl_'))))),lcl_masterBatchVars,'sorted');
end
