function [nuTBT, TCTBT, eTBT]= TBTDTparameter(r, p, n1, nc, ctn, Et, psi, l)
%% New elastic tunneling model

nt_tn1 = p.Et_res;

if Et <1e-3
    Et = 1e-3;
end

Et = Et * p.q;
DEc = 0.5 * p.q; %[J]

[a_alpha, K, K1, K2, K3, K4] = setFieldParameters(r, p, psi, 0, n1, ctn); % 1us program, Vch = -1.7V 

Fox = -K/p.e_ox/(r.r_0 + 0.5*r.t_tox1);
r_tn = r.r_tox1+(1-l)*r.t_tn1;
Fn= -(K-p.q*n1*(r.r_tox1^2-(r_tn)^2)/2)/p.e_SiON/(r_tn);

% TrapVolume=1e-15*1e-18; % [m3]
% TrapVolume = 1.75e-10^3;
TrapVolume = 4.64492e-12*1e-18;
W = 2*pi/p.h_bar*Et^2*TrapVolume;
V_c = K/p.e_ox*log(r.r_tox1/r.r_0) ...
    + K1/p.e_SiON*log((r.r_tn1-l*r.t_tn1)/r.r_tox1) ...
    + p.q*nt_tn1/4/p.e_SiON*((r.r_tox1+l*r.t_tn1)^2-r.r_tox1^2);
Etrapheight = p.PHI_B_SiON-V_c;
% if(Etrapheight < 0) 
% Etrapheight
% end
Etrapheight = max(0,p.q*Etrapheight-Et);

D = sqrt(2*p.m_e^3*.40)/pi^2/p.h_bar^3*sqrt(Etrapheight);

TC1 = exp(-4*sqrt(2*p.m_e*.5)/(3*p.h_bar*p.q)...
        * (max(0,(Et+DEc-Fn*r.t_tn1*l*p.q))^1.5-max(0,(Et+DEc-Fn*r.t_tn1*l*p.q-Fox*r.t_tox1*p.q))^1.5)/(Fox*100))...
    * exp(-4*sqrt(2*p.m_e*.4)/(3*p.h_bar*p.q)...
        * ((Et)^1.5-max(0,(Et-Fn*r.t_tn1*l*p.q))^1.5)/(Fn*100));

f = 1/(exp(-Etrapheight/p.kT)+1);
eTBT = D*W*f*TC1;
TCTBT=TC1;
nuTBT = D*W;