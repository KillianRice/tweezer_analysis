function y = asymmetric_gaussian(x, mu, sigma_left, sigma_right, offset, Amp1, Amp2)
    % x        : Input values
    % mu       : Mean (center of the Gaussian)
    % sigma_left : Standard deviation for the left side
    % sigma_right: Standard deviation for the right side
    
    % Compute the Gaussian on the left and right sides separately
    left_part = Amp1*exp(-0.5 * ((x - mu) / sigma_left).^2);
    right_part = Amp2*exp(-0.5 * ((x - mu) / sigma_right).^2);
    
    % Combine the two parts based on the direction of x
    y = left_part .* (x <= mu) + right_part .* (x > mu) + offset;
end
