function funcOut = plot_isThereAtom_tweezer(analyVar, indivDataset, avgDataset)


    % numTweezers = avgDataset.numTweezers;
% 
    % %% Plot all values together
    % x = avgDataset.imagevcoAtom;
    % y = avgDataset.meanIsThereAtom;
    % yerr = avgDataset.meanIsThereAtomStd;
% 
    % figure;
    % errorbar(x,y,yerr,...
    %     'LineStyle','none',...
    %     'Marker', 'o',...
    %     'MarkerSize', analyVar.markerSize,...
    %     'MarkerFaceColor', analyVar.COLORS(1,:),...
    %     'MarkerEdgeColor', 'k',...
    %     'Color', analyVar.COLORS(1,:));
    % title('ROI Above Threshold Percentage vs Cooling Time');
    % xlabel('Cooling Time');
    % ylabel('Signal Above Threshold Percantage');

n = avgDataset.numTweezers;

%% Original data
x = avgDataset.imagevcoAtom;
y = avgDataset.meanIsThereAtom;
yerr = avgDataset.meanIsThereAtomErr;

%% Decode x values
% x is stored as:
%     trueID.realX
%
% For example:
%     10.50 -> ID = 10, real x = 50
%     10.75 -> ID = 10, real x = 75
%     11.50 -> ID = 11, real x = 50

trueID = floor(x);
realX = round((x - trueID) * 100);

% Find all unique IDs
uniqueIDs = unique(trueID, 'stable');
numIDs = length(uniqueIDs);

% Generate a different color for each ID
colors = lines(numIDs);


%% Plot each tweezer separately
for i = 1:n

    % Take every nth value, starting at i
    x_i = realX(i:n:end);
    y_i = y(i:n:end);
    yerr_i = yerr(i:n:end);
    ID_i = trueID(i:n:end);

    figure;
    hold on;

    % Plot each ID separately so it gets its own color/legend entry
    for j = 1:numIDs

        % Points belonging to this ID
        idx = ID_i == uniqueIDs(j);

        if any(idx)
            errorbar(x_i(idx), y_i(idx), yerr_i(idx), ...
                'LineStyle', 'none', ...
                'Marker', 'o', ...
                'MarkerSize', analyVar.markerSize, ...
                'MarkerFaceColor', colors(j,:), ...
                'MarkerEdgeColor', 'k', ...
                'Color', colors(j,:), ...
                'DisplayName', sprintf('ID %g', uniqueIDs(j)));
        end

    end

    title(sprintf('Tweezer %d', i));
    xlabel('Cooling Time');
    ylabel('Signal Above Threshold Percentage');

    legend('Location', 'best');
    grid on;
    hold off;

end


%% Average every n consecutive values
numChunks = floor(length(y) / n);

% Trim to a whole number of chunks
y_trim = y(1:numChunks*n);
yerr_trim = yerr(1:numChunks*n);

% Decode the x values corresponding to the chunks
realX_trim = realX(1:numChunks*n);
ID_trim = trueID(1:numChunks*n);

% Reshape into columns:
% Each column = one chunk of n measurements
y_chunks = reshape(y_trim, n, numChunks);
yerr_chunks = reshape(yerr_trim, n, numChunks);

% The first x value of each chunk
x_avg = realX_trim(1:n:numChunks*n);

% The ID associated with each chunk
ID_avg = ID_trim(1:n:numChunks*n);

%% Calculate averaged y values
y_avg = mean(y_chunks, 1);

%% Calculate error of the averaged values
%
% Each original yerr is already the standard deviation associated
% with its corresponding averaged measurement.
%
% For the mean of n independent measurements:
%
%       sigma_mean = sqrt(sigma1^2 + ... + sigman^2) / n
%
yerr_avg = sqrt(sum(yerr_chunks.^2, 1)) / sqrt(n);


%% Plot final averaged data
figure;
hold on;

% Plot each ID separately
for j = 1:numIDs

    idx = ID_avg == uniqueIDs(j);

    if any(idx)
        errorbar(x_avg(idx), y_avg(idx), yerr_avg(idx), ...
            'LineStyle', 'none', ...
            'Marker', 'o', ...
            'MarkerSize', analyVar.markerSize, ...
            'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', 'k', ...
            'Color', colors(j,:), ...
            'DisplayName', sprintf('ID %g', uniqueIDs(j)));
    end

end

title('Average Signal Above Threshold Percentage');
xlabel('Cooling Time');
ylabel('Signal Above Threshold Percentage');

legend('Location', 'best');
grid on;
hold off;





    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
end