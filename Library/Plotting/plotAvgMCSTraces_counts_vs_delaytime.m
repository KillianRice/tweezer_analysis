function [ output_args ] = plotAvgMCSTraces_counts_vs_delaytime( analyVar,indivDataset )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[AvgSpectra, AvgSpectra_error]= deal(cell(1,indivDataset{1}.CounterAtom));
figure(analyVar.timevectorAtom(1)*10+2);
for bIndex = 1:indivDataset{1}.CounterAtom
    AvgSpectra{bIndex} = analyVar.AvgSpectra{bIndex};
    AvgSpectra_error{bIndex} = analyVar.AvgSpectra_error{bIndex};
    for DensityIndex = 1:indivDataset{1}.numDensityGroups{bIndex}
            subplot(...
                indivDataset{1}.numDensityGroups{1},...
                indivDataset{1}.CounterAtom,...
                bIndex+(DensityIndex-1)*indivDataset{1}.CounterAtom);                
%                 imagesc(delay_spectra{bIndex}(:,:,DensityIndex))
            hold on
            PeakBin = analyVar.PeakSignalBin- analyVar.roiStart+1;
            for bin = PeakBin-2:PeakBin-1
                [xticks, ~, xSpacing] = AxisTicksEngineeringForm(indivDataset{1}.timedelayOrderMatrix{bIndex});
                [yticks, ~, ySpacing] = AxisTicksEngineeringForm(AvgSpectra{bIndex}(bin,:,DensityIndex));
                y_error = AvgSpectra_error{bIndex}(bin,:,DensityIndex);
%                 plot(xticks,yticks)
                plothan=errorbar(...
                    xticks,...
                    yticks,...
                    y_error,...
                    analyVar.MARKERS2{bin},...
                    'Color',analyVar.COLORS(bin,:),...
                    'MarkerSize',analyVar.markerSize/4);
            end
            hold off
            axis square tight

            ax = gca; % current axes
            set(ax,'FontSize',8);
            set(ax,'TickDir','out');
            set(ax,'TickLength', [0.02 0.02]);
            xlabel('Field Delay (us)')
            ylabel('MCS Counts')
            title(...
                sprintf(...
                    '%s %s %s %s %s',...
                    'n0:',...
                    mat2str(DensityIndex),...
                    ', Synth:',...
                    mat2str(indivDataset{1}.synthFreq(bIndex)),...
                    'MHz'...
                )...
            )

    end
end
mtit(sprintf('%s %s %s %s','Date:', mat2str(analyVar.dataDirName),'Time:', mat2str(analyVar.timevectorAtom(1))))

     

end
