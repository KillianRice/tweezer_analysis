%% Skewed Gaussian Distribution




% Define parameters


mu = 82.0170414152363;        % Position (mean / center)
sigma =  0.122767440960753;   % Width (standard deviation)
amp = -396.93418381052;       % Peak amplitude (height)
alpha = 6.50024186035887;     % Skewness parameter (alpha > 0 right-skew, < 0 left-skew)
baseline =  1946.85023527896;

x = linspace(81.5, 82.5, 400); % X-axis grid


% Symmetric base Gaussian component
base_gaus = exp(-((x - mu).^2) / (2 * sigma^2));

% Skew component using cumulative distribution function (normcdf)
skew_factor = 2 * normcdf(alpha * (x - mu) / sigma);

skew_factor2 = 2 * normcdf(1 * (x - mu) / sigma);

% Combine with custom amplitude
y = amp * base_gaus .* skew_factor + baseline;

y2 = amp * base_gaus .* skew_factor2+ baseline;

% Plot the result
clf()
figure()
plot(x, y, 'LineWidth', 2);
hold on
plot(x,y2, 'LineWidth', 2);
grid on;
title('Arbitrary Skewed Gaussian');
xlabel('X'); ylabel('Y');

% Define the LaTeX mathematical formula string
latex_expr = '$$f(x) = 2 \cdot A \cdot \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right) \cdot \Phi\left(\alpha \frac{x-\mu}{\sigma}\right)$$';

% Add the text annotation box inside the plot bounds
text('Units', 'normalized', ...
     'Position', [0.05, 0.85], ...        % Coordinates from top-left (0 to 1)
     'String', latex_expr, ...
     'Interpreter', 'latex', ...          % Enables mathematical formatting
     'FontSize', 14, ...
     'BackgroundColor', 'w', ...          % White background
     'EdgeColor', 'k', ...                % Black border
     'Margin', 8);                        % Inside padding

hold off

