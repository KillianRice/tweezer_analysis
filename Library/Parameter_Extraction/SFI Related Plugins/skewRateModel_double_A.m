function rhat = skewRateModel_double_A(p, t)
% Double skew-peak model in count-rate units.
%
% p = [A1 mu1 sigma1 alpha1  A2 mu2 sigma2 alpha2]
% t = time vector
%
% Uses skewRateModel_A for each peak and sums them.

    A1     = p(1);  mu1    = p(2);  sigma1 = p(3);  alpha1 = p(4);
    A2     = p(5);  mu2    = p(6);  sigma2 = p(7);  alpha2 = p(8);

    rhat1 = skewRateModel_A([A1 mu1 sigma1 alpha1], t);
    rhat2 = skewRateModel_A([A2 mu2 sigma2 alpha2], t);

    rhat = rhat1 + rhat2;
end