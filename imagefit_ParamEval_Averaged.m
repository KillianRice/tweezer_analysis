function varargout = imagefit_ParamEval_Averaged(varargin)
% This program is designed to read several datafiles and corresponding background
% files, extract relevant parameters from the number distribution fits, and
% plots these normalized data sets on the same graph. The images of the
% data that is extracted are the averaged OD images of the Tweezer spots
% that contain the same dependent parameter from a scan.
%
% INPUTS:
%   varargin - variable input argument to allow passing of analysis
%              variables from analysis runner program. If not passed, the
%              program will call AnalysisData itself.
%              It is important to follow the input construction below for
%              varargin to retrieve variable data from other programs.
%              -- first  argument - analyVar
%              -- second argument - indivDataset
%              -- third argument  - avgDataset (Avg OD Image of the
%              Tweezers)
%
% OUTPUTS:
%   none
%
%
% NOTES:
%   06.10.2026 - Code edited similar to imagefit_ParamEval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%% Load variables and file data
if nargin == 0
    analyVar = AnalysisVariables;
    indivDataset = get_indiv_batch_data(analyVar);
    load(fullfile(analyVar.analyOutDir, analyVar.avgOutSubDir, 'avgDataset.mat'),'avgDataset');
else
    analyVar     = varargin{1}; % if arguments are passed analyVar must be first
    indivDataset = varargin{2}; % indivDataset must be second
    avgDataset = varargin{3};
end

if ~isfield(indivDataset{1}, 'OD_TotalCounts')
    indivDataset = param_ext_ODTotalCounts(analyVar, indivDataset);
end

if analyVar.plotFitEval && analyVar.numBasenamesAtom > 5
    disp('You are about to generate a lot of plots! Are you sure you want to continue?');
    s = input('Y/N >','s');
    while s ~= 'Y' && s ~= 'N' && s ~= 'y' && s ~= 'n'
        s = input('Y/N >','s');
    end
    if s == 'n' || s == 'N'
        return;
    end
end

fprintf('\nEvaluating averaged fit parameters...\n');

%% Modify indivDataset cell to contain fit coefficients and OD image
if analyVar.UseImages

    %% Plot averaged Tweezer ROI counts for each individual spot
    if analyVar.plotIndivTwzrCounts
        create_plot_IndivTwzr_AveragedCounts(analyVar, indivDataset, avgDataset, analyVar.plotRawCounts)
    end

    %% Plot Images of the Data
    if analyVar.fitODImage
        avgDataset = add_fit_avg_batch(analyVar,avgDataset);
    end

    if analyVar.plotRawImage
        create_plot_image_sets_ROI(analyVar, indivDataset)
    end

    %% Number Distribution Fit Evaluation
    if analyVar.plotFitEval
        % Plot cloud evolution
        create_plot_evol_avg(analyVar,avgDataset); 
        % Plot 2D fit, 1D cross section, and residuals for averaged data
        if analyVar.fitODImage
            [fit2DAxH, fit1DAxH, resAxH] = ...
                create_plot_fitEval_Averaged(analyVar, avgDataset);
        end

        %  % Standardize color limits across averaged scan
        % climMat = get_axes_prop_matrix(fit2DAxH, 'CLim');
        % ylimMat = get_axes_prop_matrix(fit1DAxH, 'YLim');
        % 
        % bestColorLim = [min(climMat(:,1)), max(climMat(:,2))];
        % best1DLim    = [-0.1, max(ylimMat(:,2))];
   % 
        %  set(evolAxH,'CLim', bestColorLim)  % apply color lim to evolution
        %  set(fit2DAxH(isgraphics(fit2DAxH)), 'CLim', bestColorLim);
        %  set(resAxH(isgraphics(resAxH)),     'CLim', bestColorLim);
        %  set(fit1DAxH(isgraphics(fit1DAxH)), 'YLim', best1DLim);
    
    end

    %% Plot averaged Tweezer ROI counts for each individual spot
    if analyVar.plotIndivTwzrCounts
        create_plot_IndivTwzr_AveragedCounts(analyVar, indivDataset, avgDataset, analyVar.plotRawCounts)
    end


    if analyVar.plotSize
        % Cloud Size - X & Y saved to avgDataset
        avgDataset = param_ext_CloudRadius_avg(analyVar,avgDataset);
        % Plotting radius of each batch
        create_plot_CloudRadius_avg(analyVar,avgDataset);
    end

     %% Tweezer Histogram
    if analyVar.plotHistogram
        % Plot the histogram of OD Total counts for all images in file
        create_plot_ODTotalCountsHistogram(analyVar,avgDataset,1)
    end


end

%% Higher order parameter fitting
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


%% Wrap Up
fclose('all'); % Close any file handles which may be open
if analyVar.SavePlotData == 1 % Output data flag in
    varargout{1} = v2struct(cat(1,'fieldNames',who()));
end