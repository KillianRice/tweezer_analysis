function plot_MCS_Integrated_SFI_roi_Spectrum(analyVar, indivDataset)

    %% Plot integrated MCS counts
    figure;
    hold on;
    for i = 1:analyVar.numBasenamesAtom

        plot(indivDataset{i}.imagevcoAtom, indivDataset{i}.sfiIntegral_roi1,...
            analyVar.MARKERS{2*i-1}, 'Color', analyVar.COLORS(2*i-1,:),...
            'MarkerSize', analyVar.markerSize,...
            'DisplayName', num2str(analyVar.timevectorAtom(i)));

        plot(indivDataset{i}.imagevcoAtom, indivDataset{i}.sfiIntegral_roi2,...
            analyVar.MARKERS{2*i}, 'Color', analyVar.COLORS(2*i,:),...
            'MarkerSize', analyVar.markerSize,...
            'DisplayName', num2str(analyVar.timevectorAtom(i)));
    end
    xlabel('GHz synth. [MHz]');
    ylabel('Total MCS Counts');
    legend('show', 'Location', 'best');
    hold off;


    %% Plot ROI1 ratio
    figure;
    hold on;
    for i = 1:analyVar.numBasenamesAtom

        plot(indivDataset{i}.imagevcoAtom, indivDataset{i}.sfiIntegral_roi1_ratio,...
            analyVar.MARKERS{i}, ...
            'Color', analyVar.COLORS(i,:), ...
            'MarkerSize', analyVar.markerSize,...
            'DisplayName', num2str(analyVar.timevectorAtom(i)));
    end
    xlabel('GHz synth. [MHz]');
    ylabel('ROI1 Ratio');
    legend('show', 'Location', 'best');
    hold off;

end