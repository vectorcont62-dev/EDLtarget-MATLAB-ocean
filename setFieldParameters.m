function x = setFieldParameters(r, p, V_ch, V_G, nt_tn1, nt_ctn, N)
q = p.q;

%% original method of solving manually
% a_alpha = 1/p.e_ox*log(r.r_tox1/r.r_0) + 1/p.e_SiON*log(r.r_tn1/r.r_tox1) + 1/p.e_ox*log(r.r_tox/r.r_tn1) + 1/p.e_n*log(r.r_ctn/r.r_tox) + 1/p.e_ox*log(r.r_box/r.r_ctn);
% K = (V_G-V_ch) / a_alpha - q*nt_tn1/2/a_alpha*(r.r_tn1^2-r.r_tox1^2)*...
%      (-r.r_tox1^2/(r.r_tn1^2 - r.r_tox1^2)/p.e_SiON*log(r.r_tn1/r.r_tox1) + 1/2/p.e_SiON + 1/p.e_ox*log(r.r_tox/r.r_tn1)+1/p.e_n*log(r.r_ctn/r.r_tox) + 1/p.e_ox*log(r.r_box/r.r_ctn))...
%      - q*nt_ctn/2/a_alpha*(r.r_ctn^2-r.r_tox^2)*...
%      (-r.r_tox^2/(r.r_ctn^2-r.r_tox^2)/p.e_n*log(r.r_ctn/r.r_tox) + 1/2/p.e_n + 1/p.e_ox*log(r.r_box/r.r_ctn));
% volN1 = pi*(r.r_tn1^2-r.r_tox1^2);
% CTG = (1/2/pi/p.e_SiON*log(r.r_tn1/(r.r_tox1+r.t_tn1*0.5)) ...
%     + 1/2/pi/p.e_ox*log(r.r_tox/r.r_tn1) + 1/2/pi/p.e_n*log(r.r_ctn/r.r_tox)...
%     + 1/2/pi/p.e_ox*log(r.r_box/r.r_ctn))^(-1);
% % CTG = (1/2/pi/p.e_ox*log(r.r_tox/r.r_tn1) + 1/2/pi/p.e_n*log(r.r_ctn/r.r_tox)...
% %     + 1/2/pi/p.e_ox*log(r.r_box/r.r_ctn))^(-1);
% volCTN = pi*(r.r_ctn^2-r.r_tox^2);
% CNG = (1/2/pi/p.e_n*log(r.r_ctn/(r.r_tox+0.5*r.t_ctn)) + 1/2/pi/p.e_ox*log(r.r_box/r.r_ctn))^(-1);
% % CNG = (1/2/pi/p.e_ox*log(r.r_box/r.r_ctn))^(-1);
% % K = 1/a_alpha*(V_G - q*nt_tn1*volN1/CTG - q*nt_ctn*volCTN/CNG);
% % K = K/5;
% K1 = K - q*nt_tn1/2*r.r_tox1^2;
% K2 = K + q*nt_tn1/2*(r.r_tn1^2-r.r_tox1^2);
% K3 = K + q*nt_tn1/2*(r.r_tn1^2-r.r_tox1^2) - q*nt_ctn/2*r.r_tox^2;
% K4 = K + q*nt_tn1/2*(r.r_tn1^2-r.r_tox1^2) + q*nt_ctn/2*(r.r_ctn^2-r.r_tox^2);

% dVt = q*nt_tn1/2*(r.r_tn1^2-r.r_tox1^2)*...
%      (-r.r_tox1^2/(r.r_tn1^2 - r.r_tox1^2)/p.e_SiON*log(r.r_tn1/r.r_tox1) + 1/2/p.e_SiON + 1/p.e_ox*log(r.r_tox/r.r_tn1)+1/p.e_n*log(r.r_ctn/r.r_tox) + 1/p.e_ox*log(r.r_box/r.r_ctn))...
%      + q*nt_ctn/2*(r.r_ctn^2-r.r_tox^2)*...
%      (-r.r_tox^2/(r.r_ctn^2-r.r_tox^2)/p.e_n*log(r.r_ctn/r.r_tox) + 1/2/p.e_n + 1/p.e_ox*log(r.r_box/r.r_ctn));

% dVtd = q*nt_tn1/2*(r.r_tn1^2-r.r_tox1^2)*...
%      (-r.r_tox1^2/(r.r_tn1^2 - r.r_tox1^2)/p.e_SiON*log(r.r_tn1/r.r_tox1) + 1/2/p.e_SiON + 1/p.e_ox*log(r.r_tox/r.r_tn1)+1/p.e_n*log(r.r_ctn/r.r_tox) + 1/p.e_ox*log(r.r_box/r.r_ctn));

%% New method of Solving by simulataneous Equation : x = A\b where Ax = b
A = zeros(N+4, N+4); % CTN slice with N + Ko1 Kn1 Ko2 KBOX
b = zeros(N+4,1);

e_o1 = p.e_ox;
e_n1 = p.e_SiON;
e_o2 = p.e_ox;
e_CTN = p.e_n;
e_BOX = p.e_ox;
r_si = r.r_si;
r_o1 = r.r_o1;
r_n1 = r.r_n1;
r_o2 = r.r_o2;
r_CTN = r.r_CTN;
r_BOX = r.r_box;

for i=1:(N+3)
     for j=1:(N+4) % simultaneous equation connection of ... 0 0 1 -1 0 0 ... except last row
          if (i==j)
               A(i,j)=1;
          elseif (i==(j-1))
               A(i,j)=-1;
          else 
               A(i,j)=0;
          end
     end
end
for k=1:(N+4) % simulataneous equation last row
     if (k==1) % Ko1 coeff
          A(N+4,k) = 1/e_o1*log(r_o1/r_si);
     elseif (k==2) % Kn1 coeff
          A(N+4,k) = 1/e_n1*log(r_n1/r_o1);
     elseif (k==3) % Ko2 coeff
          A(N+4,k) = 1/e_o2*log(r_o2/r_n1);
     elseif (k==N+4) % KBOX coeff
          A(N+4,k) = 1/e_BOX*log(r_BOX/r_CTN(N));
     elseif (k==4) % KCTN,1 coeff
          A(N+4,k) = 1/e_CTN*log(r_CTN(1)/r_o2);
     else % KCTN,2-N coeff
          A(N+4,k) = 1/e_CTN*log(r_CTN(k-3)/r_CTN(k-4));
     end
end

for l=1:(N+3) % b solution except last row
     
     if (l==1) % Ko1 coeff
          b(l) = 1/2*q*nt_tn1*r_o1^2;
     elseif (l==2) % Kn1 coeff
          b(l) = -1/2*q*nt_tn1*r_n1^2;
     elseif (l==3) % Ko2 coeff
          b(l) = 1/2*q*nt_ctn(1)*r_o2^2;
     elseif (l==N+3) % KBOX coeff
          b(l) = -1/2*q*nt_ctn(N)*r_CTN(N)^2;
     else % KCTN,2-N coeff
          b(l) = 1/2*q*nt_ctn(l-2)*r_CTN(l-3)^2 - 1/2*q*nt_ctn(l-3)*r_CTN(l-3)^2;
     end
end

b(N+4) = V_G-V_ch...
     - (q*nt_tn1/4/e_n1*(r_n1^2-r_o1^2)+q*nt_ctn(1)/4/e_CTN*(r_CTN(1)^2-r_o2^2));
for ii = 2:N
     b(N+4) = b(N+4)-q*nt_ctn(ii)/4/e_CTN*(r_CTN(ii)^2-r_CTN(ii-1)^2);
end

x = A\b;  % Solution for Ax = b >> x = A^-1 b

end