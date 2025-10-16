function [y_fit, tau, beta] = SE_fit(t, y)

% fitting by minimizing RMSE
y_max = y(1);
for i = 1:length(t)
    if(y(i) >= y_max*(exp(-1)));
        tau = t(i);
    end
end

ini_fit_rmse = RMSE(y,y_max*exp(-(t./tau).^1e-3));
for j = 2e-3:1e-3:1
    fit_rmse = RMSE(y,y_max*exp(-(t./tau).^j));
    if(fit_rmse < ini_fit_rmse)
        beta = j;
    end
    ini_fit_rmse = fit_rmse;
end

% beta = 0.2;

y_fit = y_max*exp(-(t./tau).^beta);