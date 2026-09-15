function funcOut = TweezerROILineIntegrals(analyVar, indivDataset, avgDataset)

    % ================================================================
    % ROI SETTINGS
    % ================================================================

    roiCenters = analyVar.roiCenters;
    roiRadius  = analyVar.roiRadius;

    nROIs = size(roiCenters, 1);
    nDatasets = numel(indivDataset);


    % ================================================================
    % COLORMAP
    % ================================================================

    % Blue -> white -> red
    n = 256;

    cmap = [ ...
        linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1); ...
        ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)' ...
    ];


    % ================================================================
    % LOOP THROUGH DATASETS
    % ================================================================

    for d = 1:nDatasets

        % Get filename for this dataset
        filename = analyVar.basenamevectorAtom{d};

        % Get last 8 characters of filename
        timestamp = filename(max(1, end-7):end);

        % Get stored individual ROIs
        subROI = indivDataset{d}.subROI;

        % Number of shots
        nImages = size(subROI, 1);


        % ============================================================
        % INITIALIZE LINE INTEGRAL STORAGE
        % ============================================================

        xLineIntegral = cell(nImages, nROIs);
        yLineIntegral = cell(nImages, nROIs);


        % ============================================================
        % LOOP THROUGH TWEEZERS / ROIs
        % ============================================================

        for k = 1:nROIs

            % --------------------------------------------------------
            % Determine subplot layout
            % --------------------------------------------------------

            nCols = ceil(sqrt(nImages));
            nRows = ceil(nImages / nCols);


            % ========================================================
            % FIGURE 1: INDIVIDUAL ROI IMAGES
            % ========================================================

            figure('Name', ...
                sprintf('%s - Tweezer ROI %d', timestamp, k), ...
                'Color', 'w');

            tiledlayout(nRows, nCols, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');


            % --------------------------------------------------------
            % Find common color scale for ALL shots
            % --------------------------------------------------------

            maxVal = 0;

            for i = 1:nImages

                img = subROI{i,k};

                if ~isempty(img)
                    maxVal = max(maxVal, max(abs(img(:))));
                end

            end

            lim = maxVal;


            % --------------------------------------------------------
            % Plot each shot
            % --------------------------------------------------------

            for i = 1:nImages

                nexttile;

                img = subROI{i,k};

                if isempty(img)

                    axis off;
                    title(sprintf('Shot %d - Empty', i));

                    continue;

                end

                imagesc(img, [-lim lim]);

                axis image off;

                title(sprintf('Shot %d: %g', ...
                    i, indivDataset{d}.imagevcoAtom(i)), ...
                    'FontSize', 10);

                colorbar;


                % ====================================================
                % CALCULATE LINE INTEGRALS
                % ====================================================

                % Sum over Y -> X profile
                xLineIntegral{i,k} = sum(img, 1);

                % Sum over X -> Y profile
                yLineIntegral{i,k} = sum(img, 2);

            end


            % Apply colormap
            colormap(cmap);


            % --------------------------------------------------------
            % Figure title
            % --------------------------------------------------------

            xc = roiCenters(k,1);
            yc = roiCenters(k,2);

            sgtitle(sprintf( ...
                '%s - Tweezer ROI %d  (Center: %.1f, %.1f, Radius: %.1f px)', ...
                timestamp, k, xc, yc, roiRadius));


            % ========================================================
            % FIGURE 2: X LINE INTEGRALS
            % ========================================================

            figure('Name', ...
                sprintf('%s - Tweezer ROI %d - X Integral', ...
                timestamp, k), ...
                'Color', 'w');

            tiledlayout(nRows, nCols, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');


            for i = 1:nImages

                nexttile;

                img = subROI{i,k};

                if isempty(img)

                    axis off;
                    title(sprintf('Shot %d - Empty', i));

                    continue;

                end

                % Get stored X line integral
                xIntegral = xLineIntegral{i,k};

                % X pixel coordinate
                x = 1:length(xIntegral);

                plot(x, xIntegral, ...
                    'b-', ...
                    'LineWidth', 1.5);

                grid on;

                xlabel('X (pixels)');
                ylabel('Counts');

                title(sprintf('Shot %d: %g', ...
                    i, indivDataset{d}.imagevcoAtom(i)), ...
                    'FontSize', 10);

            end

            sgtitle(sprintf( ...
                '%s - Tweezer ROI %d - X Line Integral', ...
                timestamp, k));


            % ========================================================
            % FIGURE 3: Y LINE INTEGRALS
            % ========================================================

            figure('Name', ...
                sprintf('%s - Tweezer ROI %d - Y Integral', ...
                timestamp, k), ...
                'Color', 'w');

            tiledlayout(nRows, nCols, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');


            for i = 1:nImages

                nexttile;

                img = subROI{i,k};

                if isempty(img)

                    axis off;
                    title(sprintf('Shot %d - Empty', i));

                    continue;

                end

                % Get stored Y line integral
                yIntegral = yLineIntegral{i,k};

                % Y pixel coordinate
                y = 1:length(yIntegral);

                plot(y, yIntegral, ...
                    'r-', ...
                    'LineWidth', 1.5);

                grid on;

                xlabel('Y (pixels)');
                ylabel('Counts');

                title(sprintf('Shot %d: %g', ...
                    i, indivDataset{d}.imagevcoAtom(i)), ...
                    'FontSize', 10);

            end

            sgtitle(sprintf( ...
                '%s - Tweezer ROI %d - Y Line Integral', ...
                timestamp, k));

        end


        % ============================================================
        % STORE LINE INTEGRALS IN indivDataset
        % ============================================================

        indivDataset{d}.xLineIntegral = xLineIntegral;
        indivDataset{d}.yLineIntegral = yLineIntegral;

    end


    % ================================================================
    % SAVE UPDATED indivDataset
    % ================================================================

    % Create out folder if it does not exist
    if ~exist('out', 'dir')
        mkdir('out');
    end

    % Save updated indivDataset
    save('out/indivDataset_current.mat', ...
        'indivDataset', '-v7.3');

    % Terminal output
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('indivDataset saved locally\n');
    fprintf('File: out/indivDataset_current.mat\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('New variables stored in indivDataset:\n');
    fprintf('  - subROI           : Individual ROIs for each tweezer/shot\n');
    fprintf('  - xLineIntegral    : X line integrals for each tweezer/shot\n');
    fprintf('  - yLineIntegral    : Y line integrals for each tweezer/shot\n');
    fprintf('============================================================\n');
    fprintf('\n');


    % ================================================================
    % RETURN DATA
    % ================================================================

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end
