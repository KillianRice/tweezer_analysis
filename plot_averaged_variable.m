function funcOut = plot_averaged_variable(analyVar, indivDataset, avgDataset)
%%% plot_averaged_variable.m - Brandon Torres 2025.07.29
    %%% Make a plot of the average of any quantity in indivDataset grouped
    %%% by the flags in the master batch file. To be used when doing dummy
    %%% variable scans and want to avg some quantity being measured between
    %%% different scans
    
    indVarField = 'imagevcoAtom'; % The Field of an IndivDataset that is to be plotted on the X axis
    %depVarField = 'winTotNum';
    depVarField = 'x_com';
    %depVarField = 'sfiIntegral'; % The field of an indivdataset that is to be plotted on the y axis

    %%Grab x and y values from the files and ScanIDs
    [xdata_clean, ydata_clean] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    scanIDs = analyVar.uniqScanList;

    %%Allocating space for x and y lists for plots
    x = cell(length(scanIDs));
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
            difference(num) = AvgSig(id) - AvgSig(id+1);
            wholeIDs(num) = scanIDs(id);
            diff_err(num) = sqrt(AvgSig_err(id)^2 + AvgSig_err(id+1)^2);
            num = num + 1;
        end
    end
    disp(difference)
    %%%You need to state how many unique Id's there are in the files you
    %%%are averaging since the Cell Array->avgDataset makes a structure for each
    %%%unique Id that then you index for y, yerr and x.
    %%
    %%If no unique ID's then just set UniqueID = 1
    UniqueID = 1;%%%%%
    for index = 1:UniqueID%%%%%
        avgDataset{index}.(depVarField) = y;
        avgDataset{index}.(strcat(depVarField,'_unc')) = yerr;
        avgDataset{index}.(strcat(depVarField,'_x')) = x;
    
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
        title('COM vs ID');
        xlabel('ID');
        ylabel('COM');
        hold off
    end%%%%%
    UniqueID = 1;%%%%%
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

 
    
    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
end