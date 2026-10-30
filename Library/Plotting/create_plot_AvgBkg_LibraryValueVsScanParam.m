function create_plot_AvgBkg_LibraryValueVsScanParam(analyVar, library)

xValues = library.scanParameter;
num = length(library.scanParameter);

yValues = zeros(num,1);
yErr = zeros(num,1);

for i = 1:num
    yValues(i) = mean(library.avgImageAtom);
    yErr(i) = mean(library.stdAvgImageAtom);
end

defaultAvgDataPlot(xValues,yValues,Yerr,1,analyVar)
end







function h = defaultAvgDataPlot(x,y,yerr,i,analyVar)
    h = errorbar(x,y,yerr,...
        'LineStyle','none',...
        'Marker', 'o',...
        'MarkerSize', analyVar.markerSize,...
        'MarkerFaceColor', analyVar.COLORS(i,:),...
        'MarkerEdgeColor', 'k',...
        'Color', analyVar.COLORS(i,:));
end