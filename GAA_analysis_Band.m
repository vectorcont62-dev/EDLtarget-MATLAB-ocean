function [Ec, Ev, Efield] = GAA_analysis_Band(r, p, Vch, Vg, nt_tn1, nt_ctn)

[a_alpha, K, K1, K2, K3, K4] = setFieldParameters(r, p, Vch, Vg, nt_tn1, nt_ctn);

V_c = zeros(length(r.r),1);
Ec = zeros(length(r.r),1);
Ev = zeros(length(r.r),1);
Efield = zeros(length(r.r),1);
q = p.q;
for i=1:length(r.r)
    if (r.r(i)<r.r_0)
        V_c(i) = 0;
        Ec(i) = V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_Si;
        Efield(i) = 0;

    elseif(r.r(i)>=r.r_0 && r.r(i)<r.r_tox1)
        V_c(i) = K/p.e_ox*log(r.r(i)/r.r_0);
        Ec(i) = p.PHI_B_SiO2-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - K/p.e_ox/r.r(i);    

    elseif(r.r(i)>=r.r_tox1 && r.r(i)<r.r_tn1)
        V_c(i) = K/p.e_ox*log(r.r_tox1/r.r_0) ...
            + K1/p.e_SiON*log(r.r(i)/r.r_tox1) + q*nt_tn1/4/p.e_SiON*(r.r(i)^2-r.r_tox1^2);
        Ec(i) = p.PHI_B_SiON-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiON;
        Efield(i) = - K1/p.e_SiON/r.r(i) - q*nt_tn1/2/p.e_SiON*r.r(i);

    elseif(r.r(i)>=r.r_tn1 && r.r(i)<r.r_tox)
        V_c(i) = K/p.e_ox*log(r.r_tox1/r.r_0) ...
            + K1/p.e_SiON*log(r.r_tn1/r.r_tox1) + q*nt_tn1/4/p.e_SiON*(r.r_tn1^2-r.r_tox1^2)...
            + K2/p.e_ox*log(r.r(i)/r.r_tn1);
        Ec(i) = p.PHI_B_SiO2-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - K2/p.e_ox/r.r(i);

    elseif (r.r(i)>=r.r_tox && r.r(i)<r.r_ctn)
        V_c(i) = K/p.e_ox*log(r.r_tox1/r.r_0) ...
            + K1/p.e_SiON*log(r.r_tn1/r.r_tox1) + q*nt_tn1/4/p.e_SiON*(r.r_tn1^2-r.r_tox1^2)...
            + K2/p.e_ox*log(r.r_tox/r.r_tn1)...
            + K3/p.e_n*log(r.r(i)/r.r_tox) + q*nt_ctn/4/p.e_n*(r.r(i)^2-r.r_tox^2);
        % V_c(i) = K/p.e_ox*log(r.r_tox1/r.r_0) + K1/p.e_n*log(r.r_tn1/r.r_tox1) - q*nt_tn1/4/p.e_n*(r.r_tn1^2-r.r_tox1^2)...
        %     + K2/p.e_ox*log(r.r_tox/r.r_tn1)...
        %     + K3/p.e_n*log(r.r(i)/r.r_tox) - q*nt_ctn/4/p.e_n*(r.r(i)^2-r.r_tox^2);
        Ec(i) = p.PHI_B_SiN-V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiN;
        Efield(i) = - K3/p.e_n/r.r(i) - p.q*nt_ctn/2/p.e_n*r.r(i);

    elseif (r.r(i)>= r.r_ctn && r.r(i)<r.r_box)
        V_c(i) = K/p.e_ox*log(r.r_tox1/r.r_0) ...
            + K1/p.e_SiON*log(r.r_tn1/r.r_tox1) + q*nt_tn1/4/p.e_SiON*(r.r_tn1^2-r.r_tox1^2)...
            + K2/p.e_ox*log(r.r_tox/r.r_tn1)...
            + K3/p.e_n*log(r.r_ctn/r.r_tox) + q*nt_ctn/4/p.e_n*(r.r_ctn^2-r.r_tox^2)...
            + K4/p.e_ox*log(r.r(i)/r.r_ctn);
        Ec(i) = p.PHI_B_SiO2 - V_c(i);
        Ev(i) = Ec(i) - p.E_bandgap_SiO2;
        Efield(i) = - K4/p.e_ox/r.r(i);

    else
        V_c(i) = Vg;
        Ec(i) = -Vg;
        Ev(i) = Ec(i) - 0;
        Efield(i) = 0;
    end

end

% V_diff = diff(V_c);
% Efield = -[V_diff(1,1),V_diff']'/r.dr;
end