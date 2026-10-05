% fig = openfig('C:\Users\LabUser\Downloads\junk\temp1.fig', 'invisible');
% 
% filename = 'extracted_plot_data.xlsx';
% 
% %% Find all errorbar objects
% dataObjs = findall(fig, 'Type', 'errorbar');
% 
% fprintf('Found %d errorbar objects.\n', length(dataObjs));
% 
% %% Sort objects according to their legend/display order
% dataObjs = flipud(dataObjs);
% 
% %% Determine maximum number of data points
% maxPoints = 0;
% 
% for i = 1:length(dataObjs)
% 
%     xData = dataObjs(i).XData(:);
%     yData = dataObjs(i).YData(:);
% 
%     n = min(length(xData), length(yData));
% 
%     maxPoints = max(maxPoints, n);
% 
% end
% 
% fprintf('Maximum number of points: %d\n', maxPoints);
% 
% %% Create output table
% output = table();
% 
% for i = 1:length(dataObjs)
% 
%     %% Get X and Y data
%     xData = dataObjs(i).XData(:);
%     yData = dataObjs(i).YData(:);
% 
%     n = min(length(xData), length(yData));
% 
%     xData = xData(1:n);
%     yData = yData(1:n);
% 
%     %% Get legend text
%     legendText = dataObjs(i).DisplayName;
% 
%     %% Extract ID after "="
%     equalsLocation = strfind(legendText, '=');
% 
%     if ~isempty(equalsLocation)
% 
%         groupID = strtrim( ...
%             legendText(equalsLocation(end)+1:end));
% 
%     else
% 
%         groupID = sprintf('Group_%d', i);
% 
%     end
% 
%     %% Display information
%     fprintf('Group %d: ID = %s, Points = %d\n', ...
%         i, groupID, n);
% 
%     %% Pad data with NaN
%     xData(end+1:maxPoints,1) = NaN;
%     yData(end+1:maxPoints,1) = NaN;
% 
%     %% Add columns to table
%     output.(sprintf('Group_%d_ID', i)) = ...
%         repmat({groupID}, maxPoints, 1);
% 
%     output.(sprintf('Group_%d_X', i)) = xData;
% 
%     output.(sprintf('Group_%d_Y', i)) = yData;
% 
% end
% 
% %% Delete old Excel file
% if isfile(filename)
%     delete(filename);
% end
% 
% %% Write everything to ONE sheet
% writetable(output, filename, 'Sheet', 'PlotData');
% 
% %% Close figure
% close(fig);
% 
% fprintf('\n========================================\n');
% fprintf('Finished!\n');
% fprintf('Excel file: %s\n', filename);
% fprintf('Sheet: PlotData\n');
% fprintf('Groups exported: %d\n', length(dataObjs));
% fprintf('========================================\n');

fig = openfig('C:\Users\LabUser\Downloads\junk\temp1.fig', 'invisible');

filename = 'extracted_plot_data.xlsx';

%% Find all errorbar objects
dataObjs = findall(fig, 'Type', 'errorbar');

fprintf('Found %d errorbar objects.\n', length(dataObjs));

%% Sort objects according to their legend/display order
dataObjs = flipud(dataObjs);

%% Determine maximum number of data points
maxPoints = 0;

for i = 1:length(dataObjs)

    xData = dataObjs(i).XData(:);
    yData = dataObjs(i).YData(:);

    n = min(length(xData), length(yData));

    maxPoints = max(maxPoints, n);

end

fprintf('Maximum number of points: %d\n', maxPoints);

%% Create output table
output = table();

for i = 1:length(dataObjs)

    %% Get X and Y data
    xData = dataObjs(i).XData(:);
    yData = dataObjs(i).YData(:);

    n = min(length(xData), length(yData));

    xData = xData(1:n);
    yData = yData(1:n);

    %% Get error data
    %
    % YPositiveDelta = error above the Y value
    % YNegativeDelta = error below the Y value
    %
    yPositiveError = dataObjs(i).YPositiveDelta(:);
    yNegativeError = dataObjs(i).YNegativeDelta(:);

    %% Make sure error arrays match number of data points
    yPositiveError = yPositiveError(1:min(length(yPositiveError), n));
    yNegativeError = yNegativeError(1:min(length(yNegativeError), n));

    %% Get legend text
    legendText = dataObjs(i).DisplayName;

    %% Extract ID after "="
    equalsLocation = strfind(legendText, '=');

    if ~isempty(equalsLocation)

        groupID = strtrim( ...
            legendText(equalsLocation(end)+1:end));

    else

        groupID = sprintf('Group_%d', i);

    end

    %% Display information
    fprintf('Group %d: ID = %s, Points = %d\n', ...
        i, groupID, n);

    %% Pad data with NaN
    xData(end+1:maxPoints,1) = NaN;
    yData(end+1:maxPoints,1) = NaN;

    yPositiveError(end+1:maxPoints,1) = NaN;
    yNegativeError(end+1:maxPoints,1) = NaN;

    %% Add columns to table
    output.(sprintf('Group_%d_ID', i)) = ...
        repmat({groupID}, maxPoints, 1);

    output.(sprintf('Group_%d_X', i)) = xData;

    output.(sprintf('Group_%d_Y', i)) = yData;

    output.(sprintf('Group_%d_Y_PositiveError', i)) = ...
        yPositiveError;

    output.(sprintf('Group_%d_Y_NegativeError', i)) = ...
        yNegativeError;

end

%% Delete old Excel file
if isfile(filename)
    delete(filename);
end

%% Write everything to ONE sheet
writetable(output, filename, 'Sheet', 'PlotData');

%% Close figure
close(fig);

fprintf('\n========================================\n');
fprintf('Finished!\n');
fprintf('Excel file: %s\n', filename);
fprintf('Sheet: PlotData\n');
fprintf('Groups exported: %d\n', length(dataObjs));
fprintf('========================================\n');

