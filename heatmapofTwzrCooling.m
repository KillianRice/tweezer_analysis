data = [
81.45	0.01	18.8
81.45	0.03	16.4
81.45	0.05	9.7
81.45	0.07	8
81.5	0.09	8
81.5	0.01	18.7
81.5	0.03	22.5
81.5	0.05	12.9
81.5	0.07	7.8
81.55	0.05	23.7
81.6	0.05	27.4
81.65	0.05	27.4
81.7	0.05	25.3
81.3	0.07	24
81.3	0.09	18
81.3	0.11	37.5
81.4	0.11	18.4
81.4	0.13	9.5
81.4	0.15	11.5
81.45	0.09	8.3
81.45	0.11	11
81.5	0.11	10
81.35	0.05	21
81.35	0.07	18.4
81.35	0.1	23.4
81.35	0.14	15
81.35	0.02	16.3
81.6	0.07	9.4
81.55	0.08	10.5
];
% No cooling point 81	0	20.2

x = data(:,1);
y = data(:,2);
z = data(:,3);

% Create interpolation grid
xq = linspace(min(x), max(x), 200);
yq = linspace(min(y), max(y), 200);
[X,Y] = meshgrid(xq,yq);

% Interpolate Z values
F = scatteredInterpolant(x,y,z,'natural','none');
Z = F(X,Y);

% Plot heat map
figure
imagesc(xq,yq,Z);
set(gca,'YDir','normal');

colormap(jet);
colorbar;

xlabel('Final Cooling Frequency (Mhz)');
ylabel('Set Power of Cooling Beam (RMOT) (V)');
title('Tweezer Tempeature vs Cooling power and Frequency');

hold on
scatter(x,y,80,z,'k','filled');
hold off


