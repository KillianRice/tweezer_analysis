%%Sample function for system of equations and fitting the parameters to data

%Data from scans
data = [1; 1; 1; 1; 1; 1];

initial_guess = [1 1 1];

params_est = lsqnonlin(@(p) residuals(p,data), initial_guess)

function r = residuals(params, data)
    
    p1 = params(1);
    p2 = params(2);
    p3 = params(3);

    f7 = p1^3 + p2 + p3;
    
    % your six nonlinear equations
    f1 = p1^2 + sin(p2) + p3;
    f2 = p1*exp(p2) + p3^2;
    f3 = log(p1 + p2^2) + p3;
    f4 = p1*p2 + cos(p3);
    f5 = f7;
    f6 = exp(p1) + p2^2 + p3;
    
    model = [f1; f2; f3; f4; f5; f6];
    
    % residuals
    r = model - data;

end