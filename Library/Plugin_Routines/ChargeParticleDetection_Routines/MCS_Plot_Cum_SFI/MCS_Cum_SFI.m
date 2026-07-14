function output = MCS_Cum_SFI(analyVar, indivDataset, avgDataset)
    
    % Plots the Cumulative SFI profile of a scan over a declared range of
    % the independent variable.
    % Choose the range of the independent variable over which you want
    % to accumulate SFI profiles and plot them. 

    UpperLim_imagevcoAtom = 800;      % Upper limit of independent variable named "imagevcoAtom"
    LowerLim_imagevcoAtom = 280;      % Lower limit of independent variable named "imagevcoAtom"

    for i = 1:analyVar.numBasenamesAtom  %Scan over all the file names selected in batch files.
        [roi_min, roi_max] = param_extract_sfi_roi(analyVar, indivDataset{i});  % Grab the min and max bin number in SFI scan
        indivDataset{i}.cumSFI = zeros(size(indivDataset{i}.mcsSpectra{1}));    % Prepare the matrix for x and y values
        indivDataset{i}.cumSFI(roi_min:roi_max,1) = indivDataset{i}.mcsSpectra{1}(roi_min:roi_max,1);       % Place the x values
        totsfi = zeros(size(indivDataset{i}.mcsSpectra{1}(roi_min:roi_max,2:end)));           % Dummy list that will add up all the y values for each corresponding x
        %% FOR loop over all the shots of the experiment.
        for j = 1:indivDataset{i}.CounterAtom
            % IF condition to add SFI only for a declared range above.
            %if indivDataset{i}.imagevcoAtom(j) < UpperLim_imagevcoAtom && indivDataset{i}.imagevcoAtom(j)> LowerLim_imagevcoAtom
                totsfi = totsfi + indivDataset{i}.mcsSpectra{j}(roi_min:roi_max,2:end);    %Loop for adding up the y values
            %end
        end
        indivDataset{i}.cumSFI(roi_min:roi_max,2:end) = totsfi;     % Set the y values of the matrix by the dummy variable
        %% Create the graph
        figure
        hold on
        for s = 1:size(totsfi(1,:),2)
            bar(indivDataset{i}.cumSFI(roi_min:roi_max,1),indivDataset{i}.cumSFI(roi_min:roi_max,s+1));
        end
        title(['Cumulative SFI Scan for file ',num2str(i),'; ',num2str(LowerLim_imagevcoAtom),'-',num2str(UpperLim_imagevcoAtom)], ' ')
        xlabel('Time (s)')
        ylabel('Total MCS Counts')

        %% plotting a vline for the predicted position of the signal. Only to guide the eye.
        line([analyVar.roi2_minimum*(0.0000001), analyVar.roi2_minimum*(0.0000001)], ylim, 'Color', '#A2142F', 'LineWidth', 2);
        shading flat
        hold off
    end
    
    %% Adding ScanIDs figure for looking at unique scans overlapped.
    scanIDs = analyVar.uniqScanList;
    x = cell(length(scanIDs));  % creates a nxn cell array
    y = cell(length(scanIDs));
    yerr = cell(length(scanIDs));
    for id = 1:length(scanIDs)
        x{id} = indivDataset{id}.mcsSpectra{1}(roi_min:roi_max,1);
        y{id} = zeros(size(x{id}));
        yerr{id} = zeros(size(x{id}));
        for i = 1:analyVar.numBasenamesAtom  %Scan over all the file names selected in batch files.
            if scanIDs(id) == analyVar.meanListVar(i)
                scanIDs(id);
                totsfi = zeros(size(indivDataset{i}.mcsSpectra{1}(roi_min:roi_max,2:end)));
                for j = 1:indivDataset{i}.CounterAtom
                    %if indivDataset{i}.imagevcoAtom(j) < UpperLim_imagevcoAtom && indivDataset{i}.imagevcoAtom(j)> LowerLim_imagevcoAtom
                        totsfi = totsfi + indivDataset{i}.mcsSpectra{j}(roi_min:roi_max,2:end);
                    %end
                end
            end
        end
        y{id,1} = totsfi;
    end
    
    figure;
    hold on;
    for id = 1:length(scanIDs)
        y{id};
        plot(x{id}, y{id}/sum(y{id}),'Color', analyVar.COLORS(id,:));
    end
    %sum(y{id})
    legend(num2str(scanIDs))
    title('MCS Spectra')
    xlabel('MCS bins')
    ylabel('Normalized MCS signal')
    %yscale log
    hold off
    %% End of ScanID code.
    avgDataset.('NormalizedMCSspectra') = y;
    avgDataset.('NormalizedMCSspectra_x') = x;
    output.analyVar = analyVar;
    output.indivDataset = indivDataset;
    output.avgDataset = avgDataset;
end