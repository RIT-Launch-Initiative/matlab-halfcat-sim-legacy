function [output_simulation] = output_simulation(t,thrust_lb,P_ch_psi,m_dot,P_delta_psi,figure_style)
output_simulation=figure('Name','Performance Plots vs. Time','WindowStyle','Docked');
set(gcf,'color','k');
subplot(2,2,1);
plot(t,thrust_lb,'color','red','LineWidth',1.5);
title('\color{white}THRUST');
xlabel('\color{white}Time [s]');ylabel('\color{white}Thrust [lbs]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,2);
plot(t,P_ch_psi,'color','yellow','LineWidth',1.5);
title('\color{white}CHAMBER_P');
xlabel('\color{white}Time [s]');ylabel('\color{white}Pressure [psi]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,3);
plot(t,m_dot,'color','green','LineWidth',1.5);
title('\color{white}Σm_F_L_O_W_R_A_T_E')
xlabel('\color{white}Time [s]');ylabel('\color{white}Mass Flowrate [kg/s]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,4)
plot(t,P_delta_psi,'color','blue','LineWidth',1.5);
title('\color{white}Δ_P [CHAMBER_P-OX_P]');
xlabel('\color{white}Time [s]');ylabel('\color{white}Pressure [psi]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
end

