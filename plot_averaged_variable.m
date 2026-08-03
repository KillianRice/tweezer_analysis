function funcOut = plot_averaged_variable(analyVar, indivDataset, avgDataset)
%%% plot_averaged_variable.m - Brandon Torres 2025.07.29
    %%% Make a plot of the average of any quantity in indivDataset grouped
    %%% by the flags in the master batch file. To be used when doing dummy
    %%% variable scans and want to avg some quantity being measured between
    %%% different scans
    
    indVarField = 'imagevcoAtom'; % The Field of an IndivDataset that is to be plotted on the X axis
    %depVarField = 'winTotNum';
    depVarField = 'OD_TotalCountsImg1Raw';
    %depVarField = 'sfiIntegral'; % The field of an indivdataset that is to be plotted on the y axis

    findDiff = 0; %bool that determines if you want to compare similar scans with On/Off cases
                  % Label On case as '#' on the ID slot and then OFF as
                  % '#.1'

        %% Cycle through all tweezer spots, if not using tweezers then just do once
    if analyVar.UseTweezer
        load(fullfile(analyVar.analyOutDir,'tweezerROI.mat'),'tweezerROI');
        numTweezers = size(tweezerROI.centersXY,1);
    else
        numTweezers = 1;
    end

    for tweezerNum = 1:numTweezers
        if analyVar.UseTweezer
            fprintf('plotting tweezers %d', tweezerNum);
        end

        %%Grab x and y values from the files and ScanIDs
        [xdata_clean, ydata_clean] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset, tweezerNum);
        scanIDs = analyVar.uniqScanList;
    
        %%Allocating space for x and y lists for plots
        x0 = cell(length(scanIDs));
        y = cell(length(scanIDs));
        yerr = cell(length(scanIDs));
%%  
        AvgSig = zeros(1,length(scanIDs));
        AvgSig_err = zeros(1,length(scanIDs));
    
        %%Loop over each unique ID
        for id = 1:length(scanIDs)
            %%Place holder for x values of file with unique ID
            x{id} = [];
            y{id} = [];
    
            %%Grab x and y values in each file that share same ID
            for basename = 1:analyVar.numBasenamesAtom
                if scanIDs(id) == analyVar.meanListVar(basename)
                   x{id} = union(x{id},xdata_clean{basename});
                end
            end
            for basename = 1:analyVar.numBasenamesAtom
                if scanIDs(id) == analyVar.meanListVar(basename)
                   y{id} = union(y{id},ydata_clean{basename});
                end
            end
        
    
            AvgSig(id) = mean(y{id}); % use this to average over the independent variable and get a scalar for each scanID
            AvgSig_err(id) = std(y{id})/sqrt(length(y{id})); % use this to std over the independent variable and get a scalar for each scanID
        end
        
        difference = zeros(1,floor(length(scanIDs)/2));
        wholeIDs = zeros(1,floor(length(scanIDs)/2));
        diff_err = zeros(1,floor(length(scanIDs)/2));
        num = 1;
            for id = 1:(length(scanIDs))
                if scanIDs(id) == floor(scanIDs(id))
                    if id<length(scanIDs)
                        difference(num) = AvgSig(id) - AvgSig(id+1);
                    end
                    wholeIDs(num) = scanIDs(id);
                    %diff_err(num) = sqrt(AvgSig_err(id)^2 + AvgSig_err(id+1)^2);
                    num = num + 1;
                end
            end
        
            disp(difference)
            %%%You need to state how many unique Id's there are in the files you
            %%%are averaging since the Cell Array->avgDataset makes a structure for each
            %%%unique Id that then you index for y, yerr and x.
            %%
            %%If no unique ID's then just set UniqueID = 1
            % UniqueID = length(analyVar.meanListVar);%%%%%
            % for index = 1:UniqueID%%%%%
                % avgDataset{index}.(depVarField) = y;
                % avgDataset{index}.(strcat(depVarField,'_unc')) = yerr;
                % avgDataset{index}.(strcat(depVarField,'_x')) = x;
        
            figure;
            hold on;
            for id = 1:length(scanIDs)
                errorbar(scanIDs(id), AvgSig(id), AvgSig_err(id),...
                    'LineStyle','-',...
                    'Marker', analyVar.MARKERS2(id),...
                    'MarkerSize', analyVar.markerSize,...
                    'MarkerFaceColor', analyVar.COLORS(id,:),...
                    'MarkerEdgeColor', 'none',...
                    'Color', analyVar.COLORS(id,:));
            end
            legend(num2str(scanIDs));
            title(sprintf('_Tweezer%03d',tweezerNum));
            xlabel('ID');
            ylabel('Total Averaged of Accumulated Counts');
            hold off


            % end%%%%%
            UniqueID = 1;%%%%%
            if findDiff
                for index = 1:UniqueID%%%%%
        
                    figure;
                    hold on;
                    for id = 1:length(wholeIDs)
                        errorbar(wholeIDs, difference, diff_err,...
                            'LineStyle','-',...
                            'Marker', analyVar.MARKERS2(id),...
                            'MarkerSize', analyVar.markerSize,...
                            'MarkerFaceColor', analyVar.COLORS(id,:),...
                            'MarkerEdgeColor', 'none',...
                            'Color', analyVar.COLORS(id,:));
                    end
                    legend(num2str(wholeIDs));
                    title('Difference vs ID');
                    xlabel('ID');
                    ylabel('COM Difference');
                    hold off
                end%%%%%
            end
            
    end




    % Example array (replace with your own)
    %ydata_clean;
    
    % Separate odd and even indices
    %odd_vals  = ydata_clean{1}(1:2:end);
    %even_vals = ydata_clean{1}(2:2:end);
    % data = ydata_clean{1}
    % blockSize = 10;
    % halfBlock = 5;
% 
    % odd_vals = {};
    % even_vals = {};
    % for k = 1:blockSize:numel(data)
    %     endIdx = min(k+halfBlock-1, numel(data)); % avoid overflow
    %     odd_vals = [odd_vals, data(k:endIdx)];
    % end
    % for k = halfBlock+1:blockSize:numel(data)
    %     endIdx = min(k+(halfBlock-1), numel(data)); % avoid overflow
    %     even_vals = [even_vals, data(k:endIdx)];
    % end
    % 
    % disp([odd_vals{:}])
    %     
    % % Means
    % mean_odd  = mean([odd_vals{:}]);
    % mean_even = mean([even_vals{:}]);
    % 
    % % Standard error of the mean (SEM)
    % sem_odd  = std([odd_vals{:}])  / sqrt(length([odd_vals{:}]));
    % sem_even = std([even_vals{:}]) / sqrt(length([even_vals{:}]));
    % 
    % % Display results
    % fprintf('Images with 413 On:  Mean = %.4f, SEM = %.4f\n', mean_odd, sem_odd);
    % fprintf('Images with 413 Off: Mean = %.4f, SEM = %.4f\n', mean_even, sem_even);
    % 
    % 
    % % Trim data so it's divisible by 20
    % n = floor(length(data) / blockSize) * blockSize;
    % data = data(1:n);
    % 
    % % Reshape into 20-row blocks
    % data_blocks = reshape(data, blockSize, []);
% 
    % %disp(data_blocks(1:halfBlock, :))
    % 
    % % Compute means
    % first_half_means = mean(data_blocks(1:halfBlock, :), 1);
    % second_half_means = mean(data_blocks(halfBlock+1:end, :), 1);
    % 
    % % Compute difference (last 10 - first 10)
    % diff_vals = second_half_means - first_half_means
    % 
    % % Mean and SEM
    % mean_diff = mean(diff_vals);
    % sem_diff  = std(diff_vals) / sqrt(length(diff_vals));
    % 
    % % Display results
    % fprintf('Images with 413 Off - Images with 413 On: Mean = %.4f, SEM = %.4f\n', mean_diff, sem_diff);
 
    
    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
end