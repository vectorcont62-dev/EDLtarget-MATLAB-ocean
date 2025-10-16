function [nu, TC, e] = nu_TC_relation(r, p, tspan, BETOX_mesh, nt_tn1, nt_ctn, Vch, lambda, Et_trap, BE_Conc)
nu = zeros(length(Et_trap(:,1)),2);
TC = zeros(length(Et_trap(:,1)),2);
e = zeros(length(Et_trap(:,1)),2);
for i = 1:length(Et_trap(:,1))
    nu(i,1) = Et_trap(i,1);
    TC(i,1) = Et_trap(i,1);
    e(i,1) = Et_trap(i,1);
    [nu(i,2), TC(i,2), e(i,2)] = TBTDTparameter(r, p, nt_tn1, 0, nt_ctn, Et_trap(i,1), Vch, 0.5);


end