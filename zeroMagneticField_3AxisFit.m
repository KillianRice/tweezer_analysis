function zeroMagneticField_3AxisFit()

    data = readmatrix('Analysis_2026.08.03/09102026_BFieldZero.xlsx','Sheet','Sheet1');
    
    xdata = data(:, 2:4);
    ydata = data(:, 1);

    deltaFreqModelFun = @(b, X)  sqrt( ...
        b(4).*(X(:,1)-b(1)).^2 + ...
        b(5).*(X(:,2)-b(2)).^2 + ...
        b(6).*(X(:,3)-b(3)).^2); 
    % deltaF = sqrt(Ax*(x-x0)^2 + Ay*(y-y0)^2 + Az*(z-z0)^2);

    initialguess = zeros(6,1);
    initialguess(1) = 0.115;
    initialguess(2) = 0.236;
    initialguess(3) = 1.661;
    initialguess(4) = 1;
    initialguess(5) = 1;
    initialguess(6) = 1;

    % lb = [0.1,0.1,0.1,0.1,0.1,0.1,0.01];
    % ub = [2, 2, 2, 2, 2, 2,10];

    lb = [-Inf -Inf -Inf 0 0 0];
    ub = [ Inf  Inf  Inf Inf Inf Inf];

[coeffs,resnorm,~,~,~,~,J] = lsqcurvefit(deltaFreqModelFun,initialguess,xdata,ydata,lb,ub);

n = length(ydata);
p = length(coeffs);

s2 = resnorm/(n-p);

covariance = s2 * inv(J'*J);

coeffs_err = sqrt(diag(covariance));

    disp('Zero-field voltages:')
    fprintf('X = %.3f V +/- %.4f \n', coeffs(1), full(coeffs_err(1,1)));
    fprintf('Y = %.6f V +/- %.4f \n', coeffs(2), full(coeffs_err(2,1)));
    fprintf('Z = %.6f V +/- %.4f \n', coeffs(3), full(coeffs_err(3,1)));

    disp('Fit parameters:')
    disp(coeffs)
end