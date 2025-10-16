clear; clc; %close all;
% Function to compute the modified Bessel function of the first kind
% Iv(nu, z) where nu is the order and z is the argument.

% Define the order nu and the value of z
% nu =  0;  % Example order
% z = 5;   % Example argument
% 
% % Compute the modified Bessel function Iv(nu, z)
% Iv_value = besseli(nu, z);
% 
% % Display the result
% fprintf('The value of the modified Bessel function Iv(%d, %.2f) is %.6f\n', nu, z, Iv_value);
% 
% xmax = 5; xmin = -5;
% dx = (xmax - xmin)/1e3;
% x = (xmin:dx:xmax)';
% 
% Iv_graph = besseli(nu, x);
% plot(x,Iv_graph); hold on;

%% Makram-Ebeid-Lannoo
Wt = .8*1.6e-19; Wopt = 2.5*1.6e-19; Wph = 0.2*1.6e-19;
S = (Wopt-Wt)/Wph
Fn = 1e8; Fox = 1.7e8;
% kT = 1.38e-23*298;
T = (98:10:698)';
kT = 1.38e-23*T;
nmax = ceil(Wt/Wph)-1;
n = -nmax:1:nmax; n = n';
I_prob = zeros(1,length(n));
PhononAssisted = zeros(length(n),length(kT));
TransCoeff = zeros(1,length(n));
PhononAssistedTunneling = zeros(length(kT),2);
for j = 1:length(kT)
    for i = 1:length(n)
        I_prob(i) = besseli(n(i), S/sinh(Wph/2/kT(j)));
        W = Wt + n(i)*Wph;
        TransCoeff(i) = 1.6e-19*Fn/2/sqrt(2*9.1e-31*0.5*W)*exp(-4/3/(6.6e-34/(2*pi)*1.6e-19*Fn)*sqrt(2*9.1e-31*0.5)*(W^1.5-((max(W-0.5*1e-9*Fn*1.6e-19,0)^1.5))))...
            *exp(-4/3/(6.6e-34/(2*pi)*1.6e-19*Fox)*sqrt(2*9.1e-31*0.5)*((W+0.5*1.6e-19)^1.5 - (W+0.5*1.6e-19 - 1.6e-19*Fox*2e-9)^1.5));
        PhononAssisted(j, i) = exp(Wph/2/kT(j)*n(i) - S*coth(Wph/2/kT(j))) * I_prob(i) ;
    end
    PhononAssistedTunneling(j,1) = sum(PhononAssisted(j,:));
end
%% Nasyrov-Gritsenko
% Wt = .6*1.6e-19; Wopt = 1 + Wt*1.6e-19; Wph = 0.06*1.6e-19;
% S = (Wopt-Wt)/Wph;
% Fn = 3e7; Fox = 5e7;
% m = 9.1e-31; q = 1.6e-19; hbar = 6.6e-34/(2*pi);
% % kT = 1.38e-23*298;
% T = (98:10:698)';
% kT = 1.38e-23*T;
% PhononAssisted = zeros(1,length(kT));
% for i = 1:length(kT)
%     A = 2*q*6e23*pi^0.5*hbar*Wt/(m*0.5*3e-9^2*sqrt(2*kT(i)*(Wopt-Wt)));
%     Dis_NG = exp(-(Wopt-Wt)/kT(i))*sinh(q*Fn*3e-9/2/kT(i));
%     TC = exp(-2*3e-9*sqrt(2*m*0.5*Wt)/hbar);
%     PhononAssisted(i) = A*Dis_NG*TC;
% end
%% SCLC, Mark model: Characteristic trap temperature
% e0 = 8.85e-12; q = 1.6e-19; Fn = 1e8; Fox = 1.3e8; 
% me = 9.1e-31; k = 1.38e-23; h = 6.6e-34; hbar = h/2/pi;
% Nc = (2*pi*me*0.4*T./h^2).^1.5;
% Tt = T*2; l = Tt./T; N0 = 6e17*1e6/k./Tt;
% mu = 1e-8; % [m^2/V/s]
% 
% SCLC = zeros(1,length(kT));
% for i = 1:length(SCLC)
%     SCLC(i) = mu*Nc(i)*q^(l(i)-1)*(3.9*e0/N0(i)*l(i)/(l(i)+1))^l(i)*((2*l(i)+1)/(l(i)+1))^(l(i)+1)*Fn^(l(i)+1)/3e-9^l(i);
% end

RecTemp = 1./T;
figure(5);
semilogy(1000*RecTemp, PhononAssistedTunneling(:,1)); hold on;
% semilogy(1000*RecTemp, PhononAssisted + PhononAssistedTunneling(:,1)'); hold on;
% PhononAssistedTunneling(:,1) - PhononAssistedTunneling(:, 2)
% semilogy(1000*RecTemp, SCLC); hold on;