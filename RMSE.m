function rmse_func = RMSE(y1, y2)

rmse_func = 0;
for i = 1:length(y1)
    rmse_func = rmse_func + (y1(i)-y2(i))^2;
end

