function funcOut = FitTweezerLineIntegrals(analyVar, indivDataset, avgDataset)

    %% ================================================================
    % SETTINGS
    % ================================================================

    showGuess = true;


    %% ================================================================
    % GAUSSIAN + FLAT BACKGROUND
    %
    % p = [Amplitude, Center, Sigma, Background]
    % ================================================================

    gaussFun = @(p,x) ...
        p(1) .* exp(-(x-p(2)).^2 ./ (2*p(3)^2)) + p(4);


    %% ================================================================
    % LOAD SAVED indivDataset
    % ================================================================

    dataFile = fullfile('out','indivDataset_current.mat');

    if ~exist(dataFile,'file')
        error('Could not find %s.',dataFile);
    end

    S = load(dataFile,'indivDataset');
    indivDataset = S.indivDataset;


    %% ================================================================
    % DATASET INFORMATION
    % ================================================================

    nDatasets = numel(indivDataset);
    nROIs = size(analyVar.roiCenters,1);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('Datasets contained in indivDataset\n');
    fprintf('============================================================\n');

    for d = 1:nDatasets

        filename = analyVar.basenamevectorAtom{d};

        timestamp = filename(max(1,end-7):end);

        fprintf('Dataset %d: %s\n',d,timestamp);

    end

    fprintf('============================================================\n\n');


    %% ================================================================
    % LOOP THROUGH DATASETS
    % ================================================================

    for d = 1:nDatasets

        filename = analyVar.basenamevectorAtom{d};

        timestamp = filename(max(1,end-7):end);

        xLI = indivDataset{d}.xLineIntegral;
        yLI = indivDataset{d}.yLineIntegral;

        nImages = size(xLI,1);


        % ============================================================
        % INITIALIZE STORAGE
        % ============================================================

        xFitParams = cell(nImages,nROIs);
        yFitParams = cell(nImages,nROIs);

        xFitCurve = cell(nImages,nROIs);
        yFitCurve = cell(nImages,nROIs);

        xGuessCurve = cell(nImages,nROIs);
        yGuessCurve = cell(nImages,nROIs);


        %% ============================================================
        % LOOP THROUGH TWEEZERS
        % ============================================================

        for k = 1:nROIs

            nCols = ceil(sqrt(nImages));
            nRows = ceil(nImages/nCols);


            %% ========================================================
            % X LINE INTEGRALS
            % ========================================================

            figure('Name', ...
                sprintf('%s - Tweezer %d - X Gaussian Fits', ...
                timestamp,k), ...
                'Color','w');

            tiledlayout(nRows,nCols, ...
                'TileSpacing','compact', ...
                'Padding','compact');


            for i = 1:nImages

                nexttile;

                yData = xLI{i,k};

                if isempty(yData)

                    axis off;
                    title(sprintf('Shot %d - Empty',i));
                    continue

                end

                % Force column vector
                yData = yData(:);

                % X coordinate
                xData = (1:numel(yData))';


                %% ----------------------------------------------------
                % INITIAL GUESS
                % -----------------------------------------------------

                B0 = min(yData);

                A0 = max(yData) - B0;

                [~,idx] = max(yData);

                x00 = xData(idx);

                sigma0 = max(numel(yData)/6,0.5);

                p0 = [A0,x00,sigma0,B0];

                xGuessCurve{i,k} = gaussFun(p0,xData);


                %% ----------------------------------------------------
                % PLOT DATA
                % -----------------------------------------------------

                plot(xData,yData,'b.','MarkerSize',8);

                hold on;


                %% ----------------------------------------------------
                % PLOT GUESS
                % -----------------------------------------------------

                if showGuess

                    plot(xData,xGuessCurve{i,k}, ...
                        'k--', ...
                        'LineWidth',0.75);

                end


                %% ----------------------------------------------------
                % FIT
                % -----------------------------------------------------

                try

                    % Amplitude >= 0
                    % Center inside ROI
                    % Sigma > 0
                    % Background unrestricted

                    lb = [0,1,0.1,-Inf];

                    ub = [Inf,numel(xData),numel(xData),Inf];


                    opts = optimoptions('lsqcurvefit', ...
                        'Display','off');


                    pFit = lsqcurvefit( ...
                        gaussFun, ...
                        p0, ...
                        xData, ...
                        yData, ...
                        lb, ...
                        ub, ...
                        opts);


                    % Store fit
                    xFitParams{i,k} = pFit;

                    xFitCurve{i,k} = gaussFun(pFit,xData);


                    % Plot fit
                    plot(xData,xFitCurve{i,k}, ...
                        'r-', ...
                        'LineWidth',1.5);

                    fitLabel = 'Fit';


                catch ME

                    xFitParams{i,k} = [];

                    xFitCurve{i,k} = [];

                    fitLabel = 'Fit failed';


                    % Print actual error
                    fprintf(['X fit failed: Dataset %d, ', ...
                             'Tweezer %d, Shot %d\n'], ...
                             d,k,i);

                    fprintf('    %s\n',ME.message);

                end


                hold off;

                grid on;

                xlabel('X (pixels)');
                ylabel('Counts');


                title(sprintf( ...
                    'Shot %d: %g - %s', ...
                    i, ...
                    indivDataset{d}.imagevcoAtom(i), ...
                    fitLabel), ...
                    'FontSize',10);


                if showGuess

                    if isempty(xFitCurve{i,k})

                        legend('Data','Guess', ...
                            'Location','best');

                    else

                        legend('Data','Guess','Fit', ...
                            'Location','best');

                    end

                else

                    if ~isempty(xFitCurve{i,k})

                        legend('Data','Fit', ...
                            'Location','best');

                    end

                end

            end


            sgtitle(sprintf( ...
                '%s - Tweezer %d - X Gaussian Fits', ...
                timestamp,k));


            %% ========================================================
            % Y LINE INTEGRALS
            % ========================================================

            figure('Name', ...
                sprintf('%s - Tweezer %d - Y Gaussian Fits', ...
                timestamp,k), ...
                'Color','w');

            tiledlayout(nRows,nCols, ...
                'TileSpacing','compact', ...
                'Padding','compact');


            for i = 1:nImages

                nexttile;

                yData = yLI{i,k};

                if isempty(yData)

                    axis off;
                    title(sprintf('Shot %d - Empty',i));
                    continue

                end

                % Force column vector
                yData = yData(:);

                % Coordinate
                xData = (1:numel(yData))';


                %% ----------------------------------------------------
                % INITIAL GUESS
                % -----------------------------------------------------

                B0 = min(yData);

                A0 = max(yData) - B0;

                [~,idx] = max(yData);

                x00 = xData(idx);

                sigma0 = max(numel(yData)/6,0.5);

                p0 = [A0,x00,sigma0,B0];

                yGuessCurve{i,k} = gaussFun(p0,xData);


                %% ----------------------------------------------------
                % PLOT DATA
                % -----------------------------------------------------

                plot(xData,yData,'b.','MarkerSize',8);

                hold on;


                %% ----------------------------------------------------
                % PLOT GUESS
                % -----------------------------------------------------

                if showGuess

                    plot(xData,yGuessCurve{i,k}, ...
                        'k--', ...
                        'LineWidth',0.75);

                end


                %% ----------------------------------------------------
                % FIT
                % -----------------------------------------------------

                try

                    lb = [0,1,0.1,-Inf];

                    ub = [Inf,numel(xData),numel(xData),Inf];


                    opts = optimoptions('lsqcurvefit', ...
                        'Display','off');


                    pFit = lsqcurvefit( ...
                        gaussFun, ...
                        p0, ...
                        xData, ...
                        yData, ...
                        lb, ...
                        ub, ...
                        opts);


                    % Store fit
                    yFitParams{i,k} = pFit;

                    yFitCurve{i,k} = gaussFun(pFit,xData);


                    % Plot fit
                    plot(xData,yFitCurve{i,k}, ...
                        'r-', ...
                        'LineWidth',1.5);

                    fitLabel = 'Fit';


                catch ME

                    yFitParams{i,k} = [];

                    yFitCurve{i,k} = [];

                    fitLabel = 'Fit failed';


                    % Print actual error
                    fprintf(['Y fit failed: Dataset %d, ', ...
                             'Tweezer %d, Shot %d\n'], ...
                             d,k,i);

                    fprintf('    %s\n',ME.message);

                end


                hold off;

                grid on;

                xlabel('Y (pixels)');
                ylabel('Counts');


                title(sprintf( ...
                    'Shot %d: %g - %s', ...
                    i, ...
                    indivDataset{d}.imagevcoAtom(i), ...
                    fitLabel), ...
                    'FontSize',10);


                if showGuess

                    if isempty(yFitCurve{i,k})

                        legend('Data','Guess', ...
                            'Location','best');

                    else

                        legend('Data','Guess','Fit', ...
                            'Location','best');

                    end

                else

                    if ~isempty(yFitCurve{i,k})

                        legend('Data','Fit', ...
                            'Location','best');

                    end

                end

            end


            sgtitle(sprintf( ...
                '%s - Tweezer %d - Y Gaussian Fits', ...
                timestamp,k));

        end


        %% ============================================================
        % STORE FIT RESULTS
        % ============================================================

        indivDataset{d}.xGaussianFitParams = xFitParams;
        indivDataset{d}.yGaussianFitParams = yFitParams;

        indivDataset{d}.xGaussianFitCurve = xFitCurve;
        indivDataset{d}.yGaussianFitCurve = yFitCurve;

        indivDataset{d}.xGaussianGuessCurve = xGuessCurve;
        indivDataset{d}.yGaussianGuessCurve = yGuessCurve;

    end


    %% ================================================================
    % SAVE UPDATED DATA
    % ================================================================

    save(dataFile,'indivDataset','-v7.3');


    %% ================================================================
    % RETURN
    % ================================================================

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end
