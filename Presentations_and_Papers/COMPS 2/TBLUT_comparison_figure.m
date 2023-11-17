
load('GOOD_comparison_re_data.mat','brushedData')

min_est = min(brushedData(:,1));
min_modis = min(brushedData(:,2));

max_est = max(brushedData(:,1));
max_modis = max(brushedData(:,2));

min_global = min([min_est,min_modis]);

max_global = min([max_est,max_modis]);

x = linspace(6,20,150);

avg_percent_uncert = 5.8;               % Percent! Computed value from MODIS data
MODIS_17_uncert = brushedData(:,2)*(avg_percent_uncert/100);

f = figure; plot(x,x,'k-','Linewidth',1)
hold on; grid on; grid minor
errorbar(brushedData(:,1),brushedData(:,2),MODIS_17_uncert,'vertical','m.','MarkerSize',15)
xlabel('Two Wavelength $r_{e}$ $(\mu m)$','Interpreter','latex')
ylabel('MODIS $r_{e}$ $(\mu m)$','Interpreter','latex')
set(f, 'Position', [0 0 500 500])
xlim([6 20])
ylim([6 20])
axis square


%% Create TAU_C comparison figure


load("GOOD_comparison_tauC_data.mat","brushedData1")

min_est = min(brushedData1(:,1));
min_modis = min(brushedData1(:,2));

max_est = max(brushedData1(:,1));
max_modis = max(brushedData1(:,2));

min_global = min([min_est,min_modis]);

max_global = min([max_est,max_modis]);

x = linspace(15,45,150);

avg_percent_uncert = 5.88;               % Percent! Computed value from MODIS data
MODIS_17_uncert = brushedData1(:,2)*(avg_percent_uncert/100);

f = figure; plot(x,x,'k-','Linewidth',1)
hold on; grid on; grid minor
errorbar(brushedData1(:,1),brushedData1(:,2),MODIS_17_uncert,'vertical','m.','MarkerSize',15)
xlabel('Two Wavelength $\tau_{c}$','Interpreter','latex')
ylabel('MODIS $\tau_{c}$','Interpreter','latex')
set(f, 'Position', [0 0 500 500])
xlim([15 45])
ylim([15 45])
axis square




