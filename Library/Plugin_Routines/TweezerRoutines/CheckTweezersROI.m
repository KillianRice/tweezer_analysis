function funcOut = CheckTweezersROI(analyVar, indivDataset, avgDataset)

    % ================================================================
    % ROI SETTINGS
    % ================================================================

    roiCenters = analyVar.roiCenters;
    roiRadius  = analyVar.roiRadius;

    imageSize = analyVar.matrixSize;

    nROIs = size(roiCenters, 1);
    nDatasets = numel(indivDataset);


    % ================================================================
    % LOOP THROUGH DATASETS
    % ================================================================

    for d = 1:nDatasets

        % Get raw images
        rawAtom = indivDataset{d}.rawAtomImages;
        rawBack = indivDataset{d}.rawBackImages;

        % Get filename for this dataset
        filename = analyVar.basenamevectorAtom{d};

        % Get last 8 characters of filename
        timestamp = filename(max(1, end-7):end);


        % Sum all images
        summedAtom = sum(rawAtom, 2);
        summedBack = sum(rawBack, 2);

        % Reshape into images
        summedAtom = reshape(summedAtom, imageSize);
        summedBack = reshape(summedBack, imageSize);

        % Background-subtracted summed image
        summedSub = summedAtom - summedBack;


        % ============================================================
        % LOOP THROUGH ROIs
        % ============================================================

        for k = 1:nROIs

            % ROI center
            xc = roiCenters(k, 1);
            yc = roiCenters(k, 2);

            % Pixel bounds of ROI
            xMin = max(1, floor(xc - roiRadius));
            xMax = min(imageSize(2), ceil(xc + roiRadius));

            yMin = max(1, floor(yc - roiRadius));
            yMax = min(imageSize(1), ceil(yc + roiRadius));


            % Extract ROI
            atomROI = summedAtom(yMin:yMax, xMin:xMax);
            backROI = summedBack(yMin:yMax, xMin:xMax);
            subROI  = summedSub(yMin:yMax, xMin:xMax);


            % ========================================================
            % PLOT ROI
            % ========================================================

            figure('Name', ...
                sprintf('%s - ROI %d', timestamp, k), ...
                'Color', 'w');

            tiledlayout(1, 3, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');


            % --------------------------------------------------------
            % RAW ATOM
            % --------------------------------------------------------

            nexttile;

            imagesc(atomROI);

            axis image;
            axis off;

            title('Summed Atom');

            colorbar;


            % --------------------------------------------------------
            % RAW BACKGROUND
            % --------------------------------------------------------

            nexttile;

            imagesc(backROI);

            axis image;
            axis off;

            title('Summed Background');

            colorbar;


            % --------------------------------------------------------
            % SUBTRACTED
            % --------------------------------------------------------

            nexttile;

            lim = max(abs(subROI(:)));

            imagesc(subROI, [-lim lim]);

            axis image;
            axis off;

            title('Summed Atom - Background');

            colorbar;


            % Blue -> white -> red
            n = 256;

            cmap = [ ...
                linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1); ...
                ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)' ...
            ];

            colormap(cmap);

            sgtitle(sprintf( ...
                '%s - ROI %d  (Center: %.1f, %.1f, Radius: %.1f px)', ...
                timestamp, k, xc, yc, roiRadius));

        end


        % ============================================================
        % STORE SUMMED BACKGROUND-SUBTRACTED ROIs
        % ============================================================

        summedSubROI = cell(nROIs, 1);

        for k = 1:nROIs

            % ROI center
            xc = roiCenters(k, 1);
            yc = roiCenters(k, 2);

            % Pixel bounds
            xMin = max(1, floor(xc - roiRadius));
            xMax = min(imageSize(2), ceil(xc + roiRadius));

            yMin = max(1, floor(yc - roiRadius));
            yMax = min(imageSize(1), ceil(yc + roiRadius));

            % Store summed background-subtracted ROI
            summedSubROI{k} = ...
                summedSub(yMin:yMax, xMin:xMax);

        end

        % Store in indivDataset
        indivDataset{d}.summedSubROI = summedSubROI;


        % ============================================================
        % STORE INDIVIDUAL BACKGROUND-SUBTRACTED ROIs
        % ============================================================

        nImages = size(rawAtom, 2);

        % Cell array:
        %   rows    = individual shots
        %   columns = ROIs
        subROI_individual = cell(nImages, nROIs);

        for i = 1:nImages

            % Reshape individual shot
            atomImg = reshape(rawAtom(:, i), imageSize);
            backImg = reshape(rawBack(:, i), imageSize);

            % Background subtraction BEFORE summing
            subImg = atomImg - backImg;

            for k = 1:nROIs

                % ROI center
                xc = roiCenters(k, 1);
                yc = roiCenters(k, 2);

                % Pixel bounds
                xMin = max(1, floor(xc - roiRadius));
                xMax = min(imageSize(2), ceil(xc + roiRadius));

                yMin = max(1, floor(yc - roiRadius));
                yMax = min(imageSize(1), ceil(yc + roiRadius));

                % Extract individual background-subtracted ROI
                subROI_individual{i,k} = ...
                    subImg(yMin:yMax, xMin:xMax);

            end

        end

        % Store individual ROIs in the dataset
        indivDataset{d}.subROI = subROI_individual;

    end


    % ================================================================
    % RETURN DATA
    % ================================================================

    funcOut.analyVar = analyVar;
    funcOut.indivDataset = indivDataset;
    funcOut.avgDataset = avgDataset;

end
