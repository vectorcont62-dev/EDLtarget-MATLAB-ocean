function [psi_surface] = psi_s(r, p, Vth0)
t_ox = r.t_tox1 + p.e_ox/p.e_SiON*r.t_tn1 + r.t_tox2 + p.e_ox/p.e_n*r.t_ctn + r.t_box;
Cox = p.e_ox/(r.r_0*log(1+t_ox/r.r_0));

dg_t_si = 2*r.t_si;
l = sqrt((4*p.e_si*dg_t_si+Cox*dg_t_si^2)/(8*Cox));

q = 1.6e-19; Nd = 1e15; Vfb = 0.18;
phi_0 = -Vth0 - Vfb + q*Nd/p.e_si*l^2;
psi_0 = phi_0 * (sinh(r.l_WL/l)-2*sinh(r.l_WL/2/l))/sinh(r.l_WL/l);
psi_surface = psi_0 + (-Vth0 - Vfb -psi_0)*dg_t_si^2/8/l^2;

end