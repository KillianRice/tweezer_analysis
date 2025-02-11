% Step 1: Read data from Excel sheet
filename = 'Ramp120VOutput.xlsx'; % Specify the path to your Excel file
sheet = 1; % Sheet number or sheet name
range = 'D1:E2500'; % Adjust to your data range

% Read the x and y values from the Excel file
data = readtable(filename, 'Sheet', sheet, 'Range', range);

% Extract x and y values from the table (assuming columns 1 and 2)
x = data{:, 1}; % Column 1 for x values
y = data{:, 2}; % Column 2 for y values

% Step 2: Load your existing figure
existingFig = openfig('Overlapped60-61SFI.fig'); % Specify the path to your existing figure file

% Step 3: Plot the original figure
ax = gca; % Get current axis from the existing figure

% Step 4: Add the second y-axis and plot data from Excel
yyaxis left; % Use the left y-axis for the original plot
% Ensure the original figure's plot is intact; no need to replot if it's already there.

yyaxis right; % Switch to the right y-axis for Excel data
plot(x, y, 'ro-', 'LineWidth', 2); % Plot the Excel data as red circles connected by lines

% Optional: Label the plot and customize further
xlabel('Time on Ramp (us)');
ylabel('MCS Counts');
yyaxis left; % Focus back to left axis for labeling
ylabel('Voltage of Ramp (V)');
title('SFI Scan for differing n values)');