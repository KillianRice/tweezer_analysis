function vargout = Tweezer_imagefit_ParamEval(varargin)
% This program is designed to read several datafiles and corresponding background
% files, extract relevant parameters from the number distribution fits, and
% plots these normalized data sets on the same graph.
%
% INPUTS:
%   varargin - variable input argument to allow passing of analysis
%              variables from analysis runner program. If not passed, the
%              program will call AnalysisData itself.
%              It is important to follow the input construction below for
%              varargin to retrieve variable data from other programs.
%              -- first  argument - analyVar
%              -- second argument - indivDataset
%
% OUTPUTS:
%   none
%
%
% NOTES:
%   12.10.13 - Changed name from imagefit_GaussianBimodalAndHistogramV2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load variables and file datacreate
close all
if nargin == 0 % If run without arguments
    analyVar     = TweezerAnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
    avgDataset = Tweezer_get_avg_data(analyVar,indivDataset);
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
    avgDataset = varargin{3};
end

if analyVar.numBasenamesAtom > 5
    disp('You are about to generate a lot of plots! Are you sure you want to continue?');
    s = input('Y/N >','s');
    while s ~= 'Y' && s ~= 'N' && s ~= 'y' && s ~= 'n'
        s = input('Y/N >','s');
    end
    if s == 'n' || s == 'N'
        return;
    end
end

%% Plugin Routine call for plotting and fitting
% Using the extracted parameters from above (i.e. number) we can call additional functions specified in
% AnalysisVariables to analyze and fit higher order parameters. For example, if the scans show a resonance
% lineshape, we can fit the feature and extract the resonance amplitude, position, and width. 
%
% Function called here can be specified in AnalysisVariables under Lineshape fitting but should be self
% contained. Any additional data needed should not be added to AnalysisVariables. These functions should
% perform all calculations relevant to the specified operation and should not pass data back to
% imagefit_ParamEval.
%
% The idea of this is that further analysis act as 'plugins' that do not change, modify, or obfuscate, the
% core functionality of the imagefit routine.
for lineFitIter = 1:length(analyVar.fitLineFunc)
    % Construct function handle from cell of function names
    funcHand = str2func(analyVar.fitLineFunc{lineFitIter});
    % Call function with default arguments
    funcOut = funcHand(analyVar,indivDataset,avgDataset);
    
    analyVar = funcOut.analyVar;
    indivDataset = funcOut.indivDataset;
    avgDataset = funcOut.avgDataset;
end

if exist('funcOut')    
    for index = 1:length(indivDataset)
        indivDataset{index} = funcOut.indivDataset{index};
    end
end

%% Wrap Up
fclose('all'); % Close any file handles which may be open
if analyVar.SavePlotData == 1 % Output data flag in
    vargout = v2struct(cat(1,'fieldNames',who()));
end
