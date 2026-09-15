function funcOut = CheckRawImages_v2(analyVar, indivDataset, avgDataset)

    nDatasets = numel(indivDataset);

    for d = 1:nDatasets

        % Get raw images
        rawAtom = indivDataset{d}.rawAtomImages;
        rawBack = indivDataset{d}.rawBackImages;

        % Check dimensions
        if ~isequal(size(rawAtom), size(rawBack))
            warning('Dataset %d: rawAtomImages and rawBackImages have different sizes.', d);
            continue;
        end

        nImages = size(rawAtom, 2);

        % Image dimensions
        imageSize = analyVar.matrixSize;

        nCols = ceil(sqrt(nImages));
        nRows = ceil(nImages / nCols);


        %% ============================================================
        %  FIGURE 1: RAW ATOM IMAGES
        % =============================================================

        figure('Name', sprintf('Raw Atom Images - Dataset %d', d), ...
               'Color', 'w');

        tiledlayout(nRows, nCols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        for i = 1:nImages

            nexttile;

            img = reshape(rawAtom(:, i), imageSize);

            low  = prctile(img(:), 1);
            high = prctile(img(:), 99.5);

            imagesc(img);

            axis image off;
            title(sprintf('Image %d', i));

            colorbar;
        end

        colormap hot;
        sgtitle(sprintf('Raw Atom Images - Dataset %d', d));


        %% ============================================================
        %  FIGURE 2: RAW BACKGROUND IMAGES
        % =============================================================

        figure('Name', sprintf('Raw Background Images - Dataset %d', d), ...
               'Color', 'w');

        tiledlayout(nRows, nCols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        for i = 1:nImages

            nexttile;

            img = reshape(rawBack(:, i), imageSize);

            low  = prctile(img(:), 1);
            high = prctile(img(:), 99.5);

            imagesc(img);

            axis image off;
            title(sprintf('Image %d', i), 'FontSize', 10);

            colorbar;
        end

        colormap hot;
        sgtitle(sprintf('Raw Background Images - Dataset %d', d));


        %% ============================================================
        %  FIGURE 3: SUBTRACTED IMAGES
        % =============================================================

        figure('Name', sprintf('Atom - Background - Dataset %d', d), ...
               'Color', 'w');

        tiledlayout(nRows, nCols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        % ADD: Initialize array for averaged subtracted image
        avgSubImg = zeros(imageSize);

        for i = 1:nImages

            nexttile;

            atomImg = reshape(rawAtom(:, i), imageSize);
            backImg = reshape(rawBack(:, i), imageSize);

            % Background subtraction
            subImg = atomImg - backImg;

            % ADD: Accumulate subtracted images
            avgSubImg = avgSubImg + subImg;

            lim = max(abs(subImg(:)));

            imagesc(subImg, [-lim lim]);
            
            % Blue -> white -> red
            n = 256;
            cmap = [ ...
                linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1); ...
                ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)' ...
            ];

            colormap(cmap);
            colorbar;

            % Symmetric color scale around zero
            %lim = max(abs(prctile(subImg(:), [0 100])));

            %imagesc(subImg, [-lim lim]);
            %imagesc(subImg);

            %axis image off;
            axis equal tight;
            title(sprintf('Image %d', i), 'FontSize', 10);

            colorbar;
        end

        % ADD: Calculate average subtracted image
        avgSubImg = avgSubImg / 1;

        % ADD: Plot average as the last subplot
        nexttile;

        % Use symmetric color scale around zero
        lim = max(abs(avgSubImg(:)));

        imagesc(avgSubImg, [-lim lim]);

        axis equal tight;
        title('Sum', 'FontSize', 10);

        colorbar;

        colormap(hot);

        sgtitle(sprintf('Atom - Background - Dataset %d', d));

    end
           %% ============================================================
    %  SELECT ROI CENTERS
    % =============================================================

    roiFig = figure('Name', ...
        sprintf('Select ROI Centers - Dataset %d', d), ...
        'Color', 'w');

    imagesc(avgSubImg, [-lim lim]);
    axis image;
    colormap(cmap);
    colorbar;

    title({'Select tweezer centers', ...
           'Click each center, then press ENTER when finished'});

    % Ask user for ROI radius
    answer = inputdlg( ...
        'Enter ROI radius (pixels):', ...
        'ROI Settings', ...
        [1 40], ...
        {'10'});

    if isempty(answer)
        close(roiFig);
        return;
    end

    roiRadius = str2double(answer{1});

    if isnan(roiRadius) || roiRadius <= 0
        error('ROI radius must be a positive number.');
    end

    % Select centers
    [x, y] = ginput;

    % Store centers as [x y]
    roiCenters = [x y];

    % Store in analyVar
    analyVar.roiCenters = roiCenters;
    analyVar.roiRadius = roiRadius;

    % Plot selected ROIs
    hold on;

    theta = linspace(0, 2*pi, 100);

    for k = 1:size(roiCenters, 1)

        xc = roiCenters(k, 1);
        yc = roiCenters(k, 2);

        plot(xc, yc, 'gx', ...
            'MarkerSize', 10, ...
            'LineWidth', 2);

        plot(xc + roiRadius*cos(theta), ...
             yc + roiRadius*sin(theta), ...
             'g-', ...
             'LineWidth', 1.5);

        text(xc, yc, sprintf('  %d', k), ...
            'Color', 'g', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');

    end

    hold off;

    fprintf('Selected %d ROIs\n', size(roiCenters, 1));
    fprintf('ROI radius = %.2f pixels\n', roiRadius);


    %% ============================================================
    %  SAVE ROI SETTINGS
    % =============================================================

    roiSettings.centers = roiCenters;
    roiSettings.radius = roiRadius;

    % Create "out" folder if it doesn't exist
    if ~exist('out', 'dir')
        mkdir('out');
    end

    % Save ROI settings
    roiFile = fullfile('out', ...
        sprintf('ROI_Settings_current.mat', d));

    save(roiFile, 'roiSettings');

    fprintf('ROI settings saved to:\n%s\n', roiFile);


    % Return original data
    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end
