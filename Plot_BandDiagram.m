function k = Plot_BandDiagram(Rad, E_c0, E_v0, E_field0)

figure(1);
subplot(1,2,1)
plot(Rad.r*1e7,E_c0, 'LineWidth', 2); hold on;
plot(Rad.r*1e7,E_v0, 'LineWidth', 2); hold on;

% plot(Rad.r*1e7,E_v0, 'LineWidth', 2); hold on;
legend("ConductionBand")%,"ValenceBand")
xlabel("Position [nm]")
ylabel("Energy [eV]")
title("GAA Energy Band")
low = (Rad.r_0 - 1e-7) * 1e7;
high = 1e7*Rad.r_top_channel;
xlim([low,high])
subplot(1,2,2)
plot(Rad.r*1e7,E_field0/1e6, 'LineWidth', 2); hold on;
legend("Electric Field")
xlabel("Position [nm]")
ylabel("Electric Field [MV/cm]")
title("GAA Electrid Field")

end