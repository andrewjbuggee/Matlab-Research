%% Compare my calculation of reflectance with the MODIS measurements


% By Andrew John Buggee
%%

function compare_my_reflectance_with_MODIS(modisR,R,modisInputs, pixels2use)

% ---- FOR BAND 1 ----
% grab the reflectance values estimated by my retrieval algorithm
refl_estimate_1 = [];

% grab the MODIS reflectance measurements
refl_modis_1 = modisR(:,1);


% ---- CREATE THE PLOT ----

f = figure; 
% Create axes
axes1 = axes('Parent',f,'Position',[-0.05 0.18 0.65 0.65]);
hold(axes1,'on');
% --- PLOT re FIRST ---
plot(x_r,x_r,'Parent',axes1,'Color','black','Linewidth',1)
hold on; grid on; grid minor
errorbar(est_R17,modis_R17,modis_R17_uncert,'Parent',axes1,'vertical','m.','MarkerSize',20)
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('My Estimate: $r_{e}$ $(\mu m)$','Interpreter','latex')
ylabel('MODIS Estimate: $r_{e}$ $(\mu m)$','Interpreter','latex')
title(['RMS: ',num2str(rms_diff_R17),' $\mu m$'],'Interpreter','latex')
axis square

% --- PLOT tau NEXT ---

% Create axes
axes2 = axes('Parent',f,'Position',[0.44 0.18 0.65 0.65]);
hold(axes2,'on');

plot(x_tau,x_tau,'Parent',axes2,'Color','k','Linewidth',1)
hold on; grid on; grid minor
errorbar(est_T17,modis_T17,modis_T17_uncert,'Parent',axes2,'vertical','m.','MarkerSize',20)
xlim([x_tau(1), x_tau(end)])
ylim([x_tau(1), x_tau(end)])
xlabel('My Estimate: $\tau_{c}$','Interpreter','latex')
ylabel('MODIS Estimate: $\tau_{c}$','Interpreter','latex')
title(['RMS: ',num2str(rms_diff_T17)],'Interpreter','latex')
axis square

% Create textbox
annotation('textbox',[0.111749185667752 0.678913738019169 0.131464712269273 0.0623003194888175],...
    'String',{'Bands 1 \& 7'},...
    'FontSize', 20,...
    'Interpreter','latex',...
    'FitBoxToText','on');

% set figure size
set(gcf, 'Position', [0 0 1100 700])




end