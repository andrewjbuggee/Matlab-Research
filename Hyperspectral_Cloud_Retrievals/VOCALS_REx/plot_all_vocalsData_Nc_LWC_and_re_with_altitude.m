function [] = plot_all_vocalsData_Nc_LWC_and_re_with_altitude(vocalsRex)


fnt_size_ylabel = 25;

%  Below this, create a subplot with effective radius versus time and
%  altitude versus time


figure; subplot(3,1,1)
semilogy(double(vocalsRex.time), vocalsRex.total_Nc, 'Color','k'); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex', 'Color','k', 'FontSize', fnt_size_ylabel)
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex', 'FontSize', fnt_size_ylabel)

% grab current axes
ax1 = gca;

s2 = subplot(3,1,2);
plot(double(vocalsRex.time), vocalsRex.re, 'Color', mySavedColors(18, 'fixed')); 
set(s2,'YColor', mySavedColors(18, 'fixed'));
grid on; grid minor; 
%xlabel('Time (Seconds since Takeoff)','Interpreter','latex')
ylim([0, 50])           % no need to look at effective radii larger than 50 microns

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex', 'FontSize', fnt_size_ylabel)

% change color of left y axis
yyaxis left
ylabel('$r_e$ $(\mu m)$','Interpreter','latex', 'Color',mySavedColors(18, 'fixed'), 'FontSize',fnt_size_ylabel)


% grab current axes
ax2 = gca;


s3 = subplot(3,1,3);
semilogy(double(vocalsRex.time), vocalsRex.lwc, 'Color', mySavedColors(15, 'fixed')); 
set(s3,'YColor', mySavedColors(15, 'fixed'));
grid on; grid minor; 
xlabel('Time (Seconds since Takeoff)','Interpreter','latex')
ylim([0, 50])           % no need to look at effective radii larger than 50 microns

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex', 'FontSize',fnt_size_ylabel)

% change color of left y axis
yyaxis left
ylabel('$LWC (g/m^3)$','Interpreter','latex', 'Color',mySavedColors(15, 'fixed'), 'FontSize', fnt_size_ylabel)


% grab current axes
ax3 = gca;

linkaxes([ax1 ax2 ax3],'x');


set(gcf, 'Position',[0 0 1275, 600])


end