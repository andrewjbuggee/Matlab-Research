

function [] = plot_MODIS_cloudTopHeight_alongVocalsPath(modis,vocalsRex)

% Do you want to step through each VOCALS-REx data point, or every n data
% points?
index_data2use = 1:10:length(vocalsRex.altitude);

% Create a zero array to hold all MODIS values
modis_cloudHeight_along_vocals = -1 + zeros(1, length(vocalsRex.altitude));


for ii = index_data2use
    
    % Find the minimum distance between the vocalsRex location and the
    % MODIS pixel locations
    [~,index_min] = min(sqrt((vocalsRex.latitude(ii) - modis.geo.lat).^2 + (vocalsRex.longitude(ii) - modis.geo.long).^2), [], 'all');
    
    modis_cloudHeight_along_vocals(ii) = modis.cloud.topHeight(index_min);

end

% Delete the -1's
modis_cloudHeight_along_vocals(modis_cloudHeight_along_vocals==-1) = [];

% set values where there is no cloud to nan. MODIS sets this to be -999
modis_cloudHeight_along_vocals(modis_cloudHeight_along_vocals==-999) = 0;

% Plot MODIS cloud top height along VOCALS-Rex flight path
figure; subplot(2,1,1)
plot(double(vocalsRex.time(index_data2use))./3600, modis_cloudHeight_along_vocals, 'Color',mySavedColors(2, 'fixed'))
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
title('Comparing MODIS Cloud Height with VOCALS-REx', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

yyaxis left
ylabel('Cloud Height  $(m)$','Interpreter','latex', 'Color',mySavedColors(2,'fixed'))

% grab current axes
ax1 = gca;


% Now plot the total number concentration, which is an indicator of a cloud
% being present
subplot(2,1,2)
semilogy(double(vocalsRex.time)./3600, vocalsRex.total_Nc, 'k'); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');

set(gcf, 'Position',[0 0 1275, 620])



end