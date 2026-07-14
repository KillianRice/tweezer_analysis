function funcOut = lorentzian_lineshape(analyVar, indivDataset, avgDataset)

    % Lorentzian fit
    form = @(coeffs, x) coeffs(1) ...
        .*(coeffs(3).^2)./((x-coeffs(2)).^2+coeffs(3).^2) + coeffs(4);

    saveVals = 0;

    % Toggle hard-coded vertical reference line
    plotHardCodedVLine = false;
    hardCodedVLineGHz = 127.37645901900000/8;
    hardCodedVLineMHz = hardCodedVLineGHz * 1000; % x-axis is MHz

    indVarField = 'imagevcoAtom';
    depVarField = 'sfiIntegral_roi2_ratio';

    function initialguess = x0(xdata, ydata)
        initialguess = [1, 1048.9, 0.4, 0.05];

        initialguess(1) = max(ydata) - min(ydata);
        initialguess(3) = (max(xdata) - min(xdata))/18;

        [~, position] = max(ydata);
        initialguess(2) = xdata(position);
    end

    options = struct(...
        'PlotIndivFits', false,...
        'PlotAll', false,...
        'PlotAllAvgs', true,...
        'XAxisLabel', 'WindFreak Synth (MHz)',...
        'CoeffNames', {{'Amp1.', '\mu_1.', '\sigma_1', 'Offset'}},...
        'CoeffUnits', {{'', 'MHz', 'MHz', ''}},...
        'AnnotateFunction', @myAnnotate,...
        'PlotInitialGuess', true);

    [xav, yav, yer, coefflist1, coefflist_err] = base_fit( ...
        analyVar, indivDataset, avgDataset, form, ...
        indVarField, depVarField, @x0, options);

    if plotHardCodedVLine
        hold on;
        xline(hardCodedVLineMHz, '--', ...
            sprintf('%.12f GHz', hardCodedVLineGHz), ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom');
    end

    format long g
    disp('Fit coefficients:')
    disp(coefflist1{1})
    disp('Fit coefficient errors:')
    disp(coefflist_err)

    if saveVals == 1
        filename = "MMWavePlotVals.xlsx";
        sheet    = "Sheet1";

        timeStr = strjoin(string(analyVar.timevectorAtom), ",");
        plugInStr = strjoin(string(analyVar.plugInVec), ",");

        row = {timeStr, plugInStr, ...
            coefflist1{1}(1), ...
            coefflist1{1}(2), ...
            coefflist1{1}(3), ...
            coefflist1{1}(4), ...
            '58 s', '58 d', ...
            indivDataset{1,1}.fileAtom{1}};

        if isfile(filename)
            C = readcell(filename, "Sheet", sheet);
            nextRow = size(C,1) + 1;
        else
            nextRow = 1;
        end

        range = "A" + nextRow;
        writecell(row, filename, "Sheet", sheet, "Range", range);
    end

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;
    funcOut.fitX = xav;
    funcOut.fitY = yav;
    funcOut.fitYErr = yer;
    funcOut.coefflist = coefflist1;
    funcOut.coefflist_err = coefflist_err;

end


function an = myAnnotate(coeffs, err, coeffNames, coeffUnits)

    dim = [.7 .5 .3 .3];

    nCoeff = min([numel(coeffs), numel(err), numel(coeffNames), numel(coeffUnits)]);

    strs = cell(nCoeff, 1);

    for i = 1:nCoeff
        strs{i} = sprintf('%s: %.10f +/- %.10f %s\n', ...
            coeffNames{i}, coeffs(i), err(i), coeffUnits{i});
    end

    an = annotation('textbox', dim, ...
        'String', strjoin(strs, ''), ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white');

end