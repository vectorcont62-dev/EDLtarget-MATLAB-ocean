function [Et, trap, Conc] = TrapSet_BETox(resEt, Eg, n, dt, mu, sigma)

Et = linspace(0, Eg, resEt);
Dt = zeros(1,length(Et));
%% Gaussian Trap Formation, n=1
% n = 1; %
% dt = [3e+18]; %
% mu = [1.0]; sigma = [0.05]; %

% n = 1; %1k   
% dt = [2.929e+16*0.2863/.05 + 1.687e+16*.5621/.05]; %50k
% mu = [1.0]; sigma = [0.05]; %50k

for i = 1:n
    for ii = 1:length(Et)
        Dt(ii) = Dt(ii) + dt(i)*exp(-((Et(ii)-mu(i))^2)/2/sigma(i)^2);
    end
end
Et = Et';
trap = Dt';

% trap(1) = 0;
% trap(2) = 0;

Conc = 0;
for j = 2:length(trap)
    Conc = Conc + (trap(j-1)+trap(j))/2 * Eg/resEt;
end

% startEtLevel = floor(.35/Eg*length(trap)+2);
% for j = startEtLevel:length(trap)
%     Conc = Conc + (trap(j-1)+trap(j))/2 * Eg/resEt;
% end