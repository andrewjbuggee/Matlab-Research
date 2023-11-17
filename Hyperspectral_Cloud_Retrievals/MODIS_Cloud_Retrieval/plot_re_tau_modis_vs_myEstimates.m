%% Plot the effective radius estimates versus the MODIS retrieval and the optical depth estimates in two subplots

function plot_re_tau_modis_vs_myEstimates(truth_estimate_table)


% ---------------------------------------------
% ----------- Plot Bands 1 and 7 --------------
% ---------------------------------------------

% ---- Extract droplet size ----
% extract the modis estimate and my calculation estimates
modis_R17 = truth_estimate_table.modisR17;

% The uncertainty is listed as a percent. Lets convert this to a value in
% microns
modis_R17_uncert = modis_R17.*(truth_estimate_table.modisR17_uncert./100);          % microns

% Grab my estimates
est_R17 = truth_estimate_table.estR17;

square_diffR17 = truth_estimate_table.squareDiffR17; % the absolute difference between my estimate and the modis estimate
rms_diff_R17 = sqrt(mean(square_diffR17));

% find the minimum and maximum values to create a y=x line

min_re_est = min(est_R17);
min_re_modis = min(modis_R17);

max_re_est = max(est_R17);
max_re_modis = max(modis_R17);

min_re_global = min([min_re_est,min_re_modis]);

max_re_global = max([max_re_est,max_re_modis]);

x_r = linspace((0.9 * min_re_global),(1.1*max_re_global),150);




% ----- Grab the optical depth estimates -----
% extract the modis estimate and my calculation estimates
modis_T17 = truth_estimate_table.modisT17;

% The uncertainty is listed as a percentage. Lets convert it to a
% reflectance
modis_T17_uncert = modis_T17.*(truth_estimate_table.modisT17_uncert./100);

% grab my estimates
est_T17 = truth_estimate_table.estT17;


square_diffT17 = truth_estimate_table.squareDiffT17; % the absolute difference between my estimate and the modis estimate
rms_diff_T17 = sqrt(mean(square_diffT17));


% find the minimum and maximum values to create a y=x line

min_tau_est = min(est_T17);
min_tau_modis = min(modis_T17);

max_tau_est = max(est_T17);
max_tau_modis = max(modis_T17);

min_tau_global = min([min_tau_est,min_tau_modis]);

max_tau_global = max([max_tau_est,max_tau_modis]);

x_tau = linspace((0.9 * min_tau_global),(1.1*max_tau_global),150);


% set font size for axes labels
axes_font_size = 30;
tick_labels_size = 30;

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

% set tick label font size
ax = gca;
ax.FontSize = tick_labels_size;
ax.TickLabelInterpreter = "latex";

xlabel('Two-Wavelength $r_{e}$ $(\mu m)$','Interpreter','latex', 'FontSize',axes_font_size)
ylabel('MODIS $r_{e}$ $(\mu m)$','Interpreter','latex', 'FontSize',axes_font_size)
title(['RMS: ',num2str(rms_diff_R17),' $\mu m$'],'Interpreter','latex', 'FontSize',axes_font_size)
axis square
box on

% --- PLOT tau NEXT ---

% Create axes
axes2 = axes('Parent',f,'Position',[0.44 0.18 0.65 0.65]);
hold(axes2,'on');

plot(x_tau,x_tau,'Parent',axes2,'Color','k','Linewidth',1)
hold on; grid on; grid minor
errorbar(est_T17,modis_T17,modis_T17_uncert,'Parent',axes2,'vertical','m.','MarkerSize',20)
xlim([x_tau(1), x_tau(end)])
ylim([x_tau(1), x_tau(end)])

% set tick label font size
ax = gca;
ax.FontSize = tick_labels_size;
ax.TickLabelInterpreter = "latex";

xlabel('Two-Wavelength $\tau_{c}$','Interpreter','latex', 'FontSize',axes_font_size)
ylabel('MODIS $\tau_{c}$','Interpreter','latex', 'FontSize',axes_font_size)
title(['RMS: ',num2str(rms_diff_T17)],'Interpreter','latex', 'FontSize',axes_font_size)
axis square
box on

% Create textbox
annotation('textbox',[0.111749185667752 0.678913738019169 0.131464712269273 0.0623003194888175],...
    'String',{'Bands 1 \& 7'},...
    'FontSize', 20,...
    'Interpreter','latex',...
    'FitBoxToText','on');


% set figure size
set(gcf, 'Position', [0 0 1200 800])

end