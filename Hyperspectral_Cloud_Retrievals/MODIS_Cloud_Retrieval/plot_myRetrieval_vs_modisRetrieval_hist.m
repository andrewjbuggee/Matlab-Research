%% Plot the ratio of my effective radius versus MODIs as a histogram, Do the same for optical depth

function plot_myRetrieval_vs_modisRetrieval_hist(truth_estimate_table)


% ---------------------------------------------
% ----------- Plot Bands 1 and 7 --------------
% ---------------------------------------------

% ---- Extract droplet size ----
% extract the modis estimate and my calculation estimates
modis_R17 = truth_estimate_table.modisR17;



% Grab my estimates
est_R17 = truth_estimate_table.estR17;




% ----- Grab the optical depth estimates -----
% extract the modis estimate and my calculation estimates
modis_T17 = truth_estimate_table.modisT17;



% grab my estimates
est_T17 = truth_estimate_table.estT17;




% set font size for axes labels
axes_font_size = 30;
tick_labels_size = 30;

f = figure; 
% Create axes
axes1 = axes('Parent',f,'Position',[-0.05 0.18 0.65 0.65]);
hold(axes1,'on');
% --- PLOT re FIRST ---
histogram(est_R17./modis_R17, 50)
hold on; grid on; grid minor
xlim([0.8, 1.2])

% set tick label font size
ax = gca;
ax.FontSize = tick_labels_size;
ax.TickLabelInterpreter = "latex";

xlabel('$r_{e}^{est}/r_e^{modis}$','Interpreter','latex', 'FontSize',axes_font_size)
ylabel('Counts','Interpreter','latex', 'FontSize',axes_font_size)
axis square
box on

% --- PLOT tau NEXT ---

% Create axes
axes2 = axes('Parent',f,'Position',[0.44 0.18 0.65 0.65]);
hold(axes2,'on');

histogram(est_T17./modis_T17, 50)
grid on; grid minor
xlim([0.8, 1.2])

% set tick label font size
ax = gca;
ax.FontSize = tick_labels_size;
ax.TickLabelInterpreter = "latex";

xlabel('$\tau_{c}^{est}/\tau_c^{modis}$','Interpreter','latex', 'FontSize',axes_font_size)
ylabel('Counts','Interpreter','latex', 'FontSize',axes_font_size)
axis square
box on

% Create textbox
annotation('textbox',[0.111749185667752 0.678913738019169 0.131464712269273 0.0623003194888175],...
    'String',['Bands 1 \& 7', newline, '\# pixels: ',num2str(length(est_T17))],...
    'FontSize', 20,...
    'Interpreter','latex',...
    'FitBoxToText','on');

% set figure size
set(gcf, 'Position', [0 0 1200 800])

end