function [] = plot_all_vocalsData_Nc_and_re_with_altitude(vocalsRex)


%  Below this, create a subplot with effective radius versus time and
%  altitude versus time


figure; subplot(2,1,1)
semilogy(double(vocalsRex.time)./3600, vocalsRex.total_Nc, 'Color','k'); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex', 'Color','k')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
plot(double(vocalsRex.time)./3600, vocalsRex.re, 'Color','blue'); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% change color of left y axis
yyaxis left
ylabel('$r_e$ $(\mu m)$','Interpreter','latex', 'Color','blue')


% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1275, 600])


end