function [nu0, TC0, t0] = Vth_optimize_Et_level(r, p, m, n1, ctn, psi, lambda, Et)

t0 = zeros(m,1); nu0 = zeros(m,1); TC0 = zeros(m,1);
for i = 1:m
    % mesh_mid = r.r_tox1 + (i-l)/m*r.t_tn1;
    l = (i-lambda)/m;
    [inst_nu,inst_TC,inst_e] = TBTDTparameter(r,p,n1,0,ctn,Et,psi,l);
    t0(i) = 1/inst_e; nu0(i) = inst_nu; TC0(i) = inst_TC;
end

end