function [Vt_ratio] = Vt_Vtratio_solve(r, p, tspan, BETOX_mesh, nt_tn1, nt_ctn, Vch, lambda, Et_trap, BE_Conc)

Vt_ratio = zeros(length(tspan),1);
Vt_ratio_dist = zeros(length(tspan),length(Et_trap(:,1)));
for i = 1:(length(Et_trap(:,1)))
    if i > 1
        instant_trap_density = p.E_bandgap_SiON / p.Et_res * (Et_trap(i,2) + Et_trap(i-1,2))/2;
    else
        instant_trap_density = p.E_bandgap_SiON / p.Et_res * (Et_trap(i,2));
    end
    if instant_trap_density < 1e10
        instant_trap_density = 0;
    end

    instant_trap_density = instant_trap_density/BETOX_mesh;
    [~,~,tau] = Vth_optimize_Et_level(r,p,BETOX_mesh,instant_trap_density,nt_ctn, Vch, lambda, (Et_trap(i,1)+Et_trap(max(i-1,1),1))/2);
    Vt_comp_inst = zeros(length(tspan),length(tau));
    for ii = 1:length(tau)
        for j = 1:length(tspan)
            % Vt_comp_inst(j,ii) = instant_trap_density/BE_Conc*(1-exp(-tspan(j)/tau(ii)));
            Vt_comp_inst(j,ii) = instant_trap_density*(exp(-tspan(j)/tau(ii)));
        end
    end
    
    for t = 1:length(tspan)
        Vt_ratio_dist(t,i) = sum(Vt_comp_inst(t,:));
    end
end

for tt = 1:length(tspan)
    Vt_ratio(tt) = sum(Vt_ratio_dist(tt,:));
end

% tau0 = Vth_optimize_Et_level(r,p,BETOX_mesh,nt_tn1,nt_ctn, -1.9, 0, 1.5);
% 
% Vt_comp = zeros(length(tspan),length(tau0));
% Vt_ratio = zeros(length(tspan),1);
% for i = 1:length(tau0)
%     for j = 1:length(tspan)
%         Vt_comp(j,i) = 1/length(tau0)*(1-exp(-tspan(j)/tau0(i)));
%     end
% end
% 
% for jj = 1:length(tspan)
%     Vt_ratio(jj) = sum(Vt_comp(jj,:));
% end


end