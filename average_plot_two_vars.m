function funcOut = average_plot_two_vars(analyVar, indivDataset, avgDataset)
    %%% Modified from average_plot.m - Joe Whalen 2017.12.15
    %%% Make a plot of the average of any two quantitie in indivDataset grouped
    %%% by the flags in the master batch file.
    
    indVarField = 'imagevcoAtom'; % The Field of an IndivDataset that is to be plotted on the X axis
    depVarField1 = 'numberAtom';
    depVarField2 = 'sfiIntegral'; % The field of an indivdataset that is to be plotted on the y axis
    

    xdata = cell(analyVar.numBasenamesAtom,1);
    ydata1 = cell(analyVar.numBasenamesAtom,1);
    ydata2 = cell(analyVar.numBasenamesAtom,1);
    
    for i = 1:analyVar.numBasenamesAtom
        
        xdata{i} = indivDataset{i}.(indVarField);
        ydata1{i} = indivDataset{i}.(depVarField1);
        ydata2{i} = indivDataset{i}.(depVarField2);
        
    end
    
    scanIDs = analyVar.uniqScanList;
    x1 = cell(length(scanIDs));
    y1 = cell(length(scanIDs));
    yerr1 = cell(length(scanIDs));
    x2 = cell(length(scanIDs));
    y2 = cell(length(scanIDs));
    yerr2 = cell(length(scanIDs));
    
    for id = 1:length(scanIDs)
        x1{id} = [];
        x2{id} = [];
        
        for basename = 1:analyVar.numBasenamesAtom
            if scanIDs(id) == analyVar.meanListVar(basename)
               x1{id} = union(x1{id},xdata{basename});
               x2{id} = union(x2{id},xdata{basename});
            end
        end

        y1{id} = zeros(size(x1{id}));
        yerr1{id} = zeros(size(x1{id}));
        y2{id} = zeros(size(x2{id}));
        yerr2{id} = zeros(size(x2{id}));
        
        tempy1 = zeros(size(analyVar.meanListVar));
        tempy2 = zeros(size(analyVar.meanListVar));
        for i = 1:length(x1{id})
            num=0;
            for basename = 1:analyVar.numBasenamesAtom
                disp(analyVar.numBasenamesAtom);
                if scanIDs(id) == analyVar.meanListVar(basename)
                    for j = 1:indivDataset{basename}.CounterAtom
                        if xdata{basename}(j) == x1{id}(i)
                            num = num + 1;
                            tempy1(num) = ydata1{basename}(j);
                            tempy2(num) = ydata2{basename}(j);
                        end
                    end
                end
            end
            y1{id}(i) = mean(tempy1(1:num));
            yerr1{id}(i) = std(tempy1(1:num))/5;
            y2{id}(i) = mean(tempy2(1:num));
            yerr2{id}(i) = std(tempy2(1:num))/5;
        end
    end
    
    avgDataset.(depVarField1) = y1;
    avgDataset.(strcat(depVarField1,'_unc')) = yerr1;
    avgDataset.(strcat(depVarField1,'_x')) = x1;
    avgDataset.(depVarField2) = y2;
    avgDataset.(strcat(depVarField2,'_unc')) = yerr2;
    avgDataset.(strcat(depVarField2,'_x')) = x2;

    figure;
    hold on;
    
    yyaxis left
    ylabel('Atom Number (arb.)')
    for id = 1:length(scanIDs)
        errorbar(x1{id}, y1{id}, yerr1{id},...
            'LineStyle','-',...
            'Marker', 'o',...
            'MarkerSize', analyVar.markerSize,...
            'MarkerFaceColor', analyVar.COLORS(id,:),...
            'MarkerEdgeColor', '#0072BD',...
            'Color', '#0072BD');
    end
    legend(num2str(scanIDs))
    yyaxis right
    ylabel('MCS counts')
    for id = 1:length(scanIDs)
        errorbar(x2{id}, y2{id}, yerr2{id},...
            'LineStyle','-',...
            'Marker', 's',...
            'MarkerSize', analyVar.markerSize,...
            'MarkerFaceColor', analyVar.COLORS(id,:),...
            'MarkerEdgeColor', '#A2142F',...
            'Color', '#A2142F');
    end
    legend(num2str(scanIDs))
    hold off
    
    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
end