function funcOut = average_plot(analyVar, indivDataset, avgDataset)
    %%% average_plot.m - Joe Whalen 2017.12.15
    %%% Edited on 2025.02.07 by Soumya Kanungo to include the function
    %%% 'getxy_filtered' which outputs clean datasets based on predefined
    %%% Global filters in AnalysisVariables.
    %%% Make a plot of the average of any quantity in indivDataset grouped
    %%% by the flags in the master batch file.
    
    indVarField = 'imagevcoAtom'; % The Field of an IndivDataset that is to be plotted on the X axis
    % depVarField = 'cloudRadX';
    %depVarField = 'numberAtom';
    %depVarField = 'winTotNum';
    %depVarField = 'cntrX';
    depVarField = 'sfiIntegral'; % The field of an indivdataset that is to be plotted on the y axis
    %depVarField = 'OD_TotalCounts';
    
    %[xdata, ydata] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    %[xdata_clean, ydata_clean] = getxy_filtered(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    [xdata_clean, ydata_clean] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    scanIDs = analyVar.uniqScanList;
    x = cell(length(scanIDs));
    y = cell(length(scanIDs));
    yerr = cell(length(scanIDs));
    % for id = 1:length(scanIDs)
    %     x{id} = [];
    % 
    %     for basename = 1:analyVar.numBasenamesAtom
    %         if scanIDs(id) == analyVar.meanListVar(basename)
    %            x{id} = union(x{id},xdata{basename});
    %         end
    %     end
    % 
    %     y{id} = zeros(size(x{id}));
    %     yerr{id} = zeros(size(x{id}));
    % 
    %     tempy = zeros(size(analyVar.meanListVar));
    %     for i = 1:length(x{id})
    %         num=0;
    %         for basename = 1:analyVar.numBasenamesAtom
    %             if scanIDs(id) == analyVar.meanListVar(basename)
    %                 for j = 1:indivDataset{basename}.CounterAtom
    %                     if xdata{basename}(j) == x{id}(i)
    %                         num = num + 1;
    %                         tempy(num) = ydata{basename}(j);
    %                     end
    %                 end
    %             end
    %         end
    %         num;
    %         y{id}(i) = mean(tempy(1:num));
    %         yerr{id}(i) = std(tempy(1:num))/sqrt(num);
    %     end
    % end
%% 
    AvgSig = zeros(size(scanIDs));
    AvgSig_err = zeros(size(scanIDs));
    for id = 1:length(scanIDs)
        x{id} = [];

        for basename = 1:analyVar.numBasenamesAtom
            if scanIDs(id) == analyVar.meanListVar(basename)
               x{id} = union(x{id},xdata_clean{basename});
            end
        end

        y{id} = zeros(size(x{id}));
        yerr{id} = zeros(size(x{id}));
        tempy = zeros(size(analyVar.meanListVar));
        for i = 1:length(x{id})
            num=0;
            for basename = 1:analyVar.numBasenamesAtom
                if scanIDs(id) == analyVar.meanListVar(basename)
                    for j = 1:length(xdata_clean{basename})
                        if xdata_clean{basename}(j) == x{id}(i)
                            num = num + 1;
                            tempy(num) = ydata_clean{basename}(j);
                        end
                    end
                end
            end
            num;
            y{id}(i) = mean(tempy(1:num));
            yerr{id}(i) = std(tempy(1:num))/sqrt(num);
        end
        AvgSig(id) = mean(y{id}); % use this to average over the independent variable and get a scalar for each scanID
        AvgSig_err(id) = std(y{id}); % use this to std over the independent variable and get a scalar for each scanID
    end
    
    %avgDataset.(depVarField) = y;
    %avgDataset.(strcat(depVarField,'_unc')) = yerr;
    %avgDataset.(strcat(depVarField,'_x')) = x;
    
    figure;
    hold on;
    for id = 1:length(scanIDs)
        errorbar(x{id}, y{id}, yerr{id},...
            'LineStyle','-',...
            'Marker', analyVar.MARKERS2(id),...
            'MarkerSize', analyVar.markerSize,...
            'MarkerFaceColor', analyVar.COLORS(id,:),...
            'MarkerEdgeColor', 'none',...
            'Color', analyVar.COLORS(id,:));
    end
    legend(num2str(scanIDs));
    title('689 Cooling Beam | 16 mW 532 Beam');
    xlabel('689 Cooling Beam AOM Frequency (MHz)');
    ylabel('Average Total Image Counts');
    hold off

    %% PLot
    % figure;
    % hold on;
    % errorbar(scanIDs*80e-3, AvgSig, AvgSig_err,...
    %    'LineStyle','none',...
    %    'Marker', analyVar.MARKERS2(1),...
    %    'MarkerSize', analyVar.markerSize,...
    %    'MarkerFaceColor', analyVar.COLORS(1,:),...
    %    'MarkerEdgeColor', 'none',...
    %    'Color', analyVar.COLORS(1,:));
    % title('SFI signal dependence vs time of experiment');
    % xlabel('Total time for MCS Data (ms)');
    % ylabel('Total MCS Counts');
    % hold off

    % figure;
    % hold on;
    % for id = 1:length(scanIDs)
    %     plot(x{id}, y{id}/trapz(x{id},y{id}),...
    %         'LineStyle','-',...
    %         'Marker', analyVar.MARKERS2(id),...
    %         'MarkerSize', analyVar.markerSize,...
    %         'MarkerFaceColor', analyVar.COLORS(id,:),...
    %         'MarkerEdgeColor', 'none',...
    %         'Color', analyVar.COLORS(id,:));
    % end
    % legend(num2str(scanIDs))
    % title('Trapping- Scanning Mot coil Off time - 60 1D2 Line')
    % xlabel('826 nm Synth with doubler ON [MHz]')
    % ylabel('Normalized by total signal MCS Counts')
    % hold off
    
    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    