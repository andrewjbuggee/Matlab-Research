%% Create Bohren and Clothiaux Figure 6.11


% By Andrew John Buggee

%%


function [] = plot_bohren_clothiaux_figure_6_11()





C_down = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

figure; 



% ------------------------------------------------------
% Load normal irradiances
% ------------------------------------------------------

load('2D_MC_02-Dec-2022_Wavelength_550_N-Photons_1000000_N-Layers_1_Tau0_48_SZA_0.mat','F_norm');

% Plot downward irradiance

% plot(F_norm.down,F_norm.binEdges(1:end-1)+diff(F_norm.binEdges)/2, 'Color',C_down,...
%     'LineStyle',':','LineWidth',3)

F_down = F_norm.down;
binEdges = F_norm.binEdges;

% let's plot every 3rd data point
index_down = 1:5:length(F_down);
index_edges = [index_down, length(F_down)];


plot(F_down(index_down),binEdges(index_edges(1:end-1))+diff(binEdges(index_edges))/2, 'Color',C_down,...
    'Marker','diamond','MarkerSize',7,'LineWidth',1)



hold on



% ------------------------------------------------------
% Load irradiances with solar zenith angle of 60 degrees
% ------------------------------------------------------



load('2D_MC_02-Dec-2022_Wavelength_550_N-Photons_100000_N-Layers_1_Tau0_48_SZA_60.mat','F_norm','inputs');

% Plot downward irradiance
% plot(F_norm.down,F_norm.binEdges(1:end-1)+diff(F_norm.binEdges)/2, 'Color',C_down,...
%     'LineStyle','--','LineWidth',3)


F_down = F_norm.down;
binEdges = F_norm.binEdges;

% let's plot every 3rd data point
index_down = 1:8:length(F_down);
index_edges = [index_down, length(F_down)];

plot(F_down(index_down),binEdges(index_edges(1:end-1))+diff(binEdges(index_edges))/2, 'Color',C_down,...
    'Marker','square','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',C_down)




% Create axes labels and grid

grid on; grid minor
set(gca,'YDir','reverse')
xlabel('$F/F_0$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex');
title({'Comparing Monte Carlo Simulations for ',...
    '$\theta_0= 0$ and $\theta_0 = 60$'}, 'Interpreter','latex')

% Set the yaxis tick labels
set(gca,'YTick', (0:8:48))
set(gca,'YTickLabel',(0:8:48))


% Set the size of the figure
set(gcf, 'Position',[0 0 900 500])


dim = [0.7 0.85 0 0];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
    ['$g$ = ', num2str(inputs.g)],...
    ['$r$ = ', num2str(inputs.mie.radius(1)), ' $\mu m$'],...
    ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 22;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';

legend('$F_{\downarrow}/F_0$, $\theta_0 = 0$','$F_{\downarrow}/F_0$, $\theta_0 = 60$',...
    'Interpreter','latex','Location','best','FontSize',20)


% Set x limits

xlim([0 1.2])





end