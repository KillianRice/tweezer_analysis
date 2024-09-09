% Script to run complete analysis of Neutral table data
% 
% Please be sure to have the correct date defined in AnalysisVariables
% prior to running.

% clear current figures
close all

% Load variables and file data
analyVar = AnalysisVariables;
if analyVar.numBasenamesAtom > 5
    disp(['About to analyze ' num2str(analyVar.numBasenamesAtom) ' scans, are you sure you want to continue?']);
    s = input('Y/N >','s');
    while s ~= 'Y' && s ~= 'N' && s ~= 'y' && s ~= 'n'
        s = input('Y/N >','s');
        disp(s)
    end
    if s == 'n' || s == 'N'
        return;
    end
end
indivDataset = get_indiv_batch_data(analyVar);
SumOfroiImageAtom = 0;
SumOfBackground= 0;

for basenameNum = 1:analyVar.numBasenamesAtom % number of files to analyze.
    for k = 1:indivDataset{basenameNum}.CounterAtom % number of experimental shots in each file
        disp(k);
        roiImageAtom = reshape(indivDataset{basenameNum}.AtomsCloud(:,k),[1,1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1)); 
        roiImageBackground = reshape(indivDataset{basenameNum}.BackgroundCloud(:,k),[1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1)); 
        SumOfroiImageAtom = SumOfroiImageAtom + roiImageAtom;
        %SumOfroiImageNotAtom = SumOfroiImageNotAtom + roiImageNotAtom;
        SumOfBackground = SumOfBackground + roiImageBackground;
    end
end

%Try to calculate optical density with the Intensity

%Subtract the atoms minus the background for plotting and prepare an array
%to plot transects of the image
atomsMinusBackground = SumOfroiImageAtom - SumOfBackground;
sizeOfImage = size(atomsMinusBackground);
xTransectArray= linspace(0, sizeOfImage(1) - 1, sizeOfImage(1));
yTransectArray= linspace(0, sizeOfImage(2) - 1, sizeOfImage(2));

figure
tiledlayout(2,3);

nexttile
imshow(SumOfroiImageAtom,[])
h = gca;
h.Visible = 'On';
colorbar
title("Atom Image")

nexttile
imshow(SumOfBackground, [])
colorbar
title("Background Image")
h = gca;
h.Visible = 'On';

nexttile
imshow(atomsMinusBackground, [])
colorbar
h = gca;
h.Visible = 'On';
%colormap("greyscale")
title("Atom - Background Image")

hold on


%plot a transect based on where the user clicks on the atoms-background
%image
[xTransect,yTransect] = getpts(gca)

xTransect = round(xTransect);
yTransect = round(yTransect);

xPntsForLine = zeros(sizeOfImage(1),1) + xTransect;
yPntsForLine = zeros(sizeOfImage(2),1) + yTransect;

plot(xTransectArray, yPntsForLine)
plot(xPntsForLine, yTransectArray)


nexttile
plot(xTransectArray,atomsMinusBackground(yTransect, :))
xlabel = ('x position')
ylabel = ('Proportional to Counts')
title("X-transect of Atom - Background at ", xTransect)

nexttile
plot(yTransectArray,atomsMinusBackground(:, xTransect))
xlabel = ('y position')
ylabel = ('Proportional to Counts')
title("Y-transect of Atom - Background at ", yTransect)


%imshow(SumOfroiImageNotAtom,[])
% % Plotting routine
% if analyVar.SavePlotData == 1
%     % Save data from output
%     PlotData = imagefit_ParamEval(analyVar,indivDataset);
% else
%     imagefit_ParamEval(analyVar,indivDataset);
%end