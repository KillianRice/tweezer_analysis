function create_plot_evol_avg(analyVar,avgDataset)
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

% %% Preallocate loop variables
% evolAxH = zeros(1,avgDataset.CounterAtom);
% 
% %% Loop through each avg batch
% for k = 1:avgDataset.CounterAtom
%     %% Preallocate loop variables
%     
%         % Processes all the image files in this batch/scan
%        figure(analyVar.figNum.atomEvol);
%             %% Plot Evolution
%             evolAxH(k) = subplot(avgDataset.SubPlotRows,avgDataset.SubPlotCols,k);
%             pcolor(avgDataset.All_OD_Image{k});
%             colorbar;
% 			shading flat;
% 			axis equal tight;
%             %%% Plot axis details
%             title(strcat(num2str(avgDataset.imagevcoAtom(ceil(k/avgDataset.numTweezers))),[' ' analyVar.xDataUnit]));
%     
% end
% hold on; grid on; axis on
% set(gcf,'Name',['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(1))]);
% mtit(['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(1))],'FontSize',16,'zoff',.05,'xoff',-.01)

%% Plot selected averaged OD images by imagevcoAtom value

numTweezers = avgDataset.numTweezers;

if numTweezers < 1 || numTweezers ~= round(numTweezers)
    error('avgDataset.numTweezers must be a positive integer.');
end

%% Determine the parameter values associated with the stored OD images
%
% Expected storage order:
%   Parameter 1: Tweezer 1, Tweezer 2, ..., Tweezer N
%   Parameter 2: Tweezer 1, Tweezer 2, ..., Tweezer N
%   etc.

numStoredImages = numel(avgDataset.avgODImages);

if mod(numStoredImages,numTweezers) ~= 0
    error(['The number of stored OD images (%d) is not divisible by ' ...
           'the number of tweezers (%d).'], ...
        numStoredImages,numTweezers);
end

numParameters = numStoredImages/numTweezers;

imagevcoValues = avgDataset.imagevcoAtom(:);

%% Handle either possible imagevcoAtom layout
%
% Layout A:
%   imagevcoAtom has one value per scanned parameter.
%
% Layout B:
%   imagevcoAtom has one repeated value for every tweezer image.

if numel(imagevcoValues) == numParameters

    parameterValues = imagevcoValues;

elseif numel(imagevcoValues) == numStoredImages

    parameterValues = imagevcoValues(1:numTweezers:end);

else
    error(['avgDataset.imagevcoAtom contains %d values, but expected ' ...
           'either %d values (one per parameter) or %d values ' ...
           '(one per stored image).'], ...
        numel(imagevcoValues),numParameters,numStoredImages);
end

%% Display available parameter values
fprintf('\nAvailable imagevcoAtom values:\n');

for parameterIndex = 1:numParameters
    fprintf('  %4d: %.12g\n', ...
        parameterIndex,parameterValues(parameterIndex));
end

fprintf('\nEnter one or more imagevcoAtom values.\n');
fprintf('Examples:\n');
fprintf('  1.25\n');
fprintf('  [1.25 1.50 1.75]\n');
fprintf('Press Enter without a value to cancel plotting.\n\n');

requestedValues = input('imagevcoAtom value(s) to display: ');

if isempty(requestedValues)
    fprintf('No images selected.\n');
    return;
end

if ~isnumeric(requestedValues) || any(~isfinite(requestedValues))
    error('The requested imagevcoAtom values must be finite numeric values.');
end

requestedValues = requestedValues(:);

%% Match requested values to recorded parameter values
selectedParameterIndices = zeros(size(requestedValues));

for requestNum = 1:numel(requestedValues)

    requestedValue = requestedValues(requestNum);

    % Find the closest recorded value
    [difference,parameterIndex] = min( ...
        abs(parameterValues-requestedValue));

    % Floating-point-safe matching tolerance
    matchTolerance = max( ...
        1e-12, ...
        1e-9*max(1,abs(requestedValue)));

    if difference > matchTolerance
        error(['Requested imagevcoAtom value %.12g was not found.\n' ...
               'The closest available value is %.12g.'], ...
            requestedValue,parameterValues(parameterIndex));
    end

    selectedParameterIndices(requestNum) = parameterIndex;
end

% Remove repeated selections while preserving entered order
selectedParameterIndices = unique( ...
    selectedParameterIndices,'stable');

%% Choose a readable subplot arrangement for the tweezers
numPlotRows = floor(sqrt(numTweezers));
numPlotRows = max(numPlotRows,1);

numPlotCols = ceil(numTweezers/numPlotRows);

while numPlotRows*numPlotCols < numTweezers
    numPlotRows = numPlotRows+1;
end

%% Determine global color limits for all selected images

globalMin = inf;
globalMax = -inf;

for selectionNum = 1:numel(selectedParameterIndices)

    parameterIndex = selectedParameterIndices(selectionNum);

    for tweezerNum = 1:numTweezers

        imageIndex = (parameterIndex-1)*numTweezers + tweezerNum;

        if imageIndex > numStoredImages
            continue;
        end

        img = avgDataset.avgODImages{imageIndex};

        globalMin = min(globalMin, min(img(:), [], 'omitnan'));
        globalMax = max(globalMax, max(img(:), [], 'omitnan'));

    end
end

%% Make one figure for each selected parameter value
for selectionNum = 1:numel(selectedParameterIndices)

    parameterIndex = selectedParameterIndices(selectionNum);
    parameterValue = parameterValues(parameterIndex);

    figH = figure( ...
        'Name',sprintf( ...
            'Cloud Evolution: imagevcoAtom = %.12g %s', ...
            parameterValue,analyVar.xDataUnit), ...
        'NumberTitle','off');

    tiledH = tiledlayout( ...
        figH,numPlotRows,numPlotCols, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    axH = gobjects(numTweezers,1);

    for tweezerNum = 1:numTweezers

        %% Convert parameter/tweezer indices to stored-image index
        imageIndex = ...
            (parameterIndex-1)*numTweezers + tweezerNum;

        if imageIndex > numStoredImages
            warning(['Missing stored image for parameter index %d, ' ...
                     'tweezer %d.'], ...
                parameterIndex,tweezerNum);
            continue;
        end

        thisODImage = avgDataset.avgODImages{imageIndex};

        axH(tweezerNum) = nexttile(tiledH,tweezerNum);

        pcolor(axH(tweezerNum),thisODImage);
        shading(axH(tweezerNum),'flat');

        axis(axH(tweezerNum),'image');
        axis(axH(tweezerNum),'tight');

        colorbar(axH(tweezerNum));
        clim(axH(tweezerNum), [globalMin globalMax]);

        xlabel(axH(tweezerNum),'Column pixel');
        ylabel(axH(tweezerNum),'Row pixel');

        title(axH(tweezerNum),sprintf( ...
            'Tweezer %d\n%s = %.12g %s', ...
            tweezerNum, ...
            'imagevcoAtom', ...
            parameterValue, ...
            analyVar.xDataUnit));

        grid(axH(tweezerNum),'on');
        box(axH(tweezerNum),'on');
    end

    title(tiledH,sprintf( ...
        ['Averaged OD Images\n' ...
         '%.12g %s'], ...
        parameterValue,analyVar.xDataUnit), ...
        'FontSize',16);
end