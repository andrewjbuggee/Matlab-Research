function [] = plot_all_vocalsData_Nc_and_re_with_altitude(vocalsRex)


%  Below this, create a subplot with effective radius versus time and
%  altitude versus time


figure; subplot(2,1,1)
semilogy(double(vocalsRex.time), vocalsRex.total_Nc, 'Color','k'); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex', 'Color','k')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

s2 = subplot(2,1,2);
plot(double(vocalsRex.time), vocalsRex.re, 'Color', mySavedColors(18, 'fixed')); 
set(s2,'YColor', mySavedColors(4, 'fixed'));
grid on; grid minor; 
xlabel('Time (Seconds since Takeoff)','Interpreter','latex')
ylim([0, 50])           % no need to look at effective radii larger than 50 microns

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% change color of left y axis
yyaxis left
ylabel('$r_e$ $(\mu m)$','Interpreter','latex', 'Color',mySavedColors(18, 'fixed'))


% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1275, 600])


end