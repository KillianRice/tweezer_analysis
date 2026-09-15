function funcOut = FittedCoeffs_plot(analyVar, indivDataset, avgDataset)
    
    coeff2Plot = 3;  %%% Write Coeff number based on fit being used
    disp(avgDataset.avg_fit_coeffs_twzr)

    xlabeltext = "B field Voltage";
    ylabeltext = "Width";
% 
    % indVarField = 'imagevcoAtom';
    % depVarField = 'sfiIntegral';

    scanIDs = avgDataset.avg_fit_coeffs_twzr.scanIDs;
    %scanIDs = analyVar.uniqScanList;
    coeffs = avgDataset.avg_fit_coeffs_twzr.coeffs;
    %[xdata_clean, ydata_clean] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset);
    %[blue,blue_err,purple,purple_err] = get_daq_averages(analyVar, indivDataset);
    %[num, num_err, tx, ~, ty, ~] = get_num_temp_averages(analyVar, indivDataset);

    %ydata_clean
    uncs = avgDataset.avg_fit_coeffs_twzr.uncs;
    numTweezers = avgDataset.avg_fit_coeffs_twzr.numTweezers;
    
    % coeffs = cellfun(@sum, ydata_clean);
    % coeffs_err = sqrt(coeffs);

    %coeffs = coeffs./blue./purple
    % coeffs = coeffs'./blue./purple./num
    % coeffs_err = coeffs_err'./blue./purple./num

    % coeffs_std = coeffs .* sqrt( ...
    % (coeffs_err ./ coeffs).^2 + ...
    % (blue_err ./ blue).^2 + ...
    % (purple_err ./ purple).^2 + ...
    % (num_err ./ num).^2 );
    for tweezerNum = 1:numTweezers
        tempCoeffs = [];
        tempUncs = [];
        for i = 1:length(scanIDs{tweezerNum})
            tempCoeffs(end+1,1) = coeffs{tweezerNum}{i}(coeff2Plot);
            tempUncs(end+1,1) = uncs{tweezerNum}{i}(coeff2Plot);
        end
    
        figure
        hold on
        errorbar(scanIDs{tweezerNum}, tempCoeffs, tempUncs,...
            'LineStyle','none',...
            'Marker', 'o',...
            'MarkerSize', analyVar.markerSize,...
            'MarkerFaceColor', analyVar.COLORS(tweezerNum,:),...
            'MarkerEdgeColor', 'k',...
            'Color', analyVar.COLORS(tweezerNum,:));
         legend(num2str(tweezerNum));
         title(sprintf('Fti Coeff %d for Tweezer ROI %d',coeff2Plot,tweezerNum));
         hold off
    
    end

    avg_coeffs = [];
    avg_uncs = [];
    for i = 1:length(scanIDs{1})
        tempCoeff = 0;
        tempUnc = 0;
        for num = 1:numTweezers
            tempCoeff = tempCoeff + coeffs{num}{i}(coeff2Plot);
            tempUnc = tempUnc + uncs{num}{i}(coeff2Plot);
        end
        avg_coeffs(end+1,1) = tempCoeff / numTweezers;
        avg_uncs(end+1,1) = tempUnc / numTweezers;
    end
    % figure
    % hold on
    % errorbar(scanIDs, coeffs, coeffs_std,...
    %     'LineStyle','none',...
    %     'Marker', 'o');
    %  xlabel(xlabeltext,'Interpreter','none');
    %  ylabel(ylabeltext,'Interpreter','none');
    %  legend('Data');
    %  %title(sprintf('Avg Fit Coeff %d for all Tweezers',coeff2Plot));
    %  title(sprintf('Signal = Counts/N/I_{461}/I_{413})'));
    %  annotation('textbox', [0.2 0.7 0.2 0.1], ...
    % 'String', 'Front plate = -1500 V \newlineCM = -205 V', ...
    % 'FitBoxToText', 'on');
    %  hold off


    figure
    hold on
    errorbar(scanIDs{1}, avg_coeffs, avg_uncs,...
        'LineStyle','none',...
        'Marker', 'o',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(tweezerNum,:),...
        'MarkerEdgeColor', 'k',...
        'Color', analyVar.COLORS(tweezerNum,:));
     xlabel(xlabeltext,'Interpreter','none');
     ylabel(ylabeltext,'Interpreter','none');
     %legend(num2str(1),'Data','Fit');
     title(sprintf('Avg Fit Coeff %d for all Tweezers',coeff2Plot));
     hold off


    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;