%% Master file for creating plots to use in COMPS 2 Report


% By Andrew John Buggee


%%  Create plot showing the change in reflectance versus optical depth. 
% Show the analytical two-stream solution as well?

clear variables

% Load the data sets!

% Changing optical depth
Tau = [2:2:50];



% Create the reflectance vector in a for loop
R = [];
% Create absorption vector
A = [];
% Create Transmission vector
T = [];
% Create the sum of all components vector
S = [];

for nn = 1:length(Tau)
    
    % grab file name

    filename{nn} = ['2D_MC_19-Feb-2023_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_',num2str(Tau(nn)),'_SZA_45.mat'];

    load(filename{nn})

    % add to reflectance vector
    R(nn) = final_state.scatter_out_top/inputs.N_photons;

    % add to reflectance vector
    A(nn) = final_state.absorbed/inputs.N_photons;

    % Add to transmit vector
    T(nn) = final_state.scatter_out_bottom/inputs.N_photons;

    % Compute the sum for each optical depth. These should all be 1!
    S(nn) = R(nn) + A(nn) + T(nn);

end


% Create figure
figure;

% Create the first x-axis plot
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1, R, Tau, 'Color',mySavedColors(8, 'fixed'))
ax1.XLabel.String = 'Reflectance';
ax1.XLabel.Color = mySavedColors(8, 'fixed');
ax1.XColor = mySavedColors(8, 'fixed');
ax1.XLabel.Interpreter = 'latex';ax1.YDir = "reverse";
ax1.XLim = [0 0.6];
grid on; grid minor
ax1.YLabel.String = 'Cloud Optical Depth';
ax1.YLabel.Color = 'k';
ax1.YLabel.Interpreter = 'latex';
ax1.YDir = "reverse";


% Create the second x axis plot

ax2 = axes(t);
plot(ax2, A,Tau,'Color',mySavedColors(9, 'fixed'))
ax2.XAxisLocation = 'top';
%ax2.YAxisLocation = 'right';

% Set the color of the axes object to 'none' so that the underlying plot is visible.
% Turn off the plot boxes to prevent the box edges from obscuring the x- and y-axes.

ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax2.YDir = "reverse";

% Turn off y tick labels
ax2.YTickLabel = [];

ax2.XLabel.String = 'Absorptance';
ax2.XLabel.Color = mySavedColors(9, 'fixed');
ax2.XColor = mySavedColors(9, 'fixed');
ax2.XLabel.Interpreter = 'latex';ax1.YDir = "reverse";
ax2.XLim = [0 0.6];


% Create textbox with simulation properties

% Textbox
dim = [0.168 0.183549713134768 0.181382743835449 0.440450286865234];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['$N_{layers}$ = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$\mu_0$ = ', num2str(round(cosd(inputs.solar_zenith_angle), 3))],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';

% Set figure size
set(gcf,'Position', [0 0 1000 500])






% ----------------------------------------------------
% OR with 1 yaxis of Optical depth and 1 common x axis

% Create figure
figure;

% Create the first x-axis plot

plot(R, Tau, 'Color',mySavedColors(8, 'fixed'))
hold on
plot(A, Tau, 'Color',mySavedColors(9, 'fixed'))
plot(T, Tau, 'Color',mySavedColors(10, 'fixed'))
set(gca,'YDir', 'reverse')

ylabel('Cloud Optical Depth', 'Interpreter','latex');
xlabel('Fraction of All Photons ($n/N_{total}$)', 'Interpreter','latex')
grid on; grid minor
legend('Reflectance','Absorptance','Transmittance','Location','best', 'Interpreter','latex',...
    'FontSize',20)


% Create textbox with simulation properties

% Textbox
dim = [0.706307100085543 0.236958094278979 0.181382743835449 0.367041905721029];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['$N_{layers}$ = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$\mu_0$ = ', num2str(round(cosd(inputs.solar_zenith_angle), 3))],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';

% Set figure size
set(gcf,'Position', [0 0 1000 600])





% -----------------------------------
% Or with optical depth on the x axis

% Create figure
figure;

% Create the first x-axis plot

plot(Tau, R, 'Color',mySavedColors(8, 'fixed'))
hold on
plot(Tau, A, 'Color',mySavedColors(9, 'fixed'))
plot(Tau, T, 'Color',mySavedColors(10, 'fixed'))

xlabel('Cloud Optical Depth', 'Interpreter','latex');
ylabel('Fraction of All Photons ($n/N_{total}$)', 'Interpreter','latex')
grid on; grid minor
legend('Reflectance','Absorptance','Transmittance','Location','bestoutside', 'Interpreter','latex',...
    'FontSize',20)


% Create textbox with simulation properties

% Textbox
dim = [0.754307100085542 0.253624760945646 0.181382743835449 0.367041905721029];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['$N_{layers}$ = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$\mu_0$ = ', num2str(round(cosd(inputs.solar_zenith_angle), 3))],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';

% Set figure size
set(gcf,'Position', [0 0 1000 600])
