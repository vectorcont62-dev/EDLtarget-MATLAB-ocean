function [tau, beta] = compact_tau(r, p, n1, nc, ctn, Et, sigma)


[~, ~, e_max] = TBTDTparameter(r, p, n1, nc, ctn, Et - 1.24*sigma, 0, 0); % 1.28: Gaussian table 10% point
[~, ~, e_min] = TBTDTparameter(r, p, n1, nc, ctn, Et + 1.24*sigma, 0, 1);



% [~, ~, e_] = TBTDTparameter(r, p, n1, nc, ctn, 2.8, 0, .99)

tau_min = 1/e_max;
tau_max = 1/e_min;
tau = sqrt(1/e_max/e_min);
% tau = (tau_min)^(exp(-1))*(tau_max)^(1-(exp(-1)));
% Delta = tau_max/tau_min
% t = 0.5;
% tau = (tau_min)^(t)*(tau_max)^(1-t);
beta = log(log(1/0.9))/log(sqrt(tau_min/tau_max));
% beta = log(log(1/0.1))/log(sqrt(tau_max/tau_min));
% beta = .5*log(log(1/0.95))/log(sqrt(tau_min/tau_max));