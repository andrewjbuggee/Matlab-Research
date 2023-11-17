%% Recreating Bohren and Clothiaux Figure 6.10 Monte Carlo Plot


% By Andrew John Buggee


%% Make Plot

function [] = plot_bohren_clothiaux_figure_6_10()

% -----------------------------------------------
% ***** LOAD RESULTS FOR OPTICAL DEPTH OF 2 *****
% -----------------------------------------------

sim_results = load('./Monte_Carlo_Simulation_Results/2D_MC_02-Dec-2022_Wavelength_550_N-Photons_100000_N-Layers_1_Tau0_2_SZA_0.mat');

% -------------------------
% ***** Unpack Inputs *****
% -------------------------

if sim_results.inputs.mie.integrate_over_size_distribution==true
    ssa = sim_results.inputs.ssa_avg;
    g = sim_results.inputs.g_avg;
else
    ssa = sim_results.inputs.ssa;
    g = sim_results.inputs.g;
end


tau_y_lower_limit = sim_results.inputs.tau_y_lower_limit;
tau_y_upper_limit = sim_results.inputs.tau_y_upper_limit;

% Define the albedo of the lower boundary in our model
albedo_maxTau = sim_results.inputs.albedo_maxTau;


binEdges = sim_results.F_norm.binEdges;



% --------------------------------------------------------
% ------- PLOT TWO STREAM THEORY F_UP AND F_DOWN ---------
% --------------------------------------------------------


% K is defined in Bohren and Clothiaux (eq. 5.70)
K = sqrt((1 - ssa)*(1 - g*ssa));
% Define the reflectivity at the top of our layer, the photons that
% scatter out the cloud top
R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

% Define the constants
A = (R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


% lets plot some range of tau
tau = linspace(tau_y_lower_limit,tau_y_upper_limit,100);

C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

% Define linewidth
LW = 3;

figure;
set(gcf, 'Position', [0 0 600 630]);
s1 = subplot(3,1,1);

plot(photon_fraction_up(tau),tau,'Color',C1,'LineStyle','-','LineWidth',LW)
hold on;
plot(photon_fraction_down(tau),tau, 'Color',C2, 'LineStyle','-','LineWidth',LW)


% -------------------------------------------------------------
% ------- PLOT MONTE CARLO SIMULATION F_UP AND F_DOWN ---------
% -------------------------------------------------------------

plot(sim_results.F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle',':','LineWidth',LW)
plot(sim_results.F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle',':','LineWidth',LW)
grid on; grid minor
set(gca,'YDir','reverse')
ylabel('Optical Depth','Interpreter','latex', 'FontSize',17);

ylim([tau_y_lower_limit, tau_y_upper_limit])
xlim([0, 1.2])

% remove x axis tick label
set(gca,'XTickLabel',[])

% Set the yaxis tick labels
set(gca,'YTick', (0:0.4:2))
set(gca,'YTickLabel',(0:0.4:2))

set(gca,'FontSize',12)

% Print the legend last so the algorithm doesn't place it in front
% of the text box
legend('$F_{\uparrow}/F_0$ 2-stream','$F_{\downarrow}/F_0$ 2-stream',...
    '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ 2D Monte Carlo','Interpreter','latex',...
    'Location','best','FontSize',18)

%%



% -------------------------------------------------
% ***** LOAD RESULTS FOR OPTICAL DEPTH OF 6.5 *****
% -------------------------------------------------

sim_results = load('./Monte_Carlo_Simulation_Results/2D_MC_02-Dec-2022_Wavelength_550_N-Photons_100000_N-Layers_1_Tau0_6.5_SZA_0.mat');

% -------------------------
% ***** Unpack Inputs *****
% -------------------------

if sim_results.inputs.mie.integrate_over_size_distribution==true
    ssa = sim_results.inputs.ssa_avg;
    g = sim_results.inputs.g_avg;
else
    ssa = sim_results.inputs.ssa;
    g = sim_results.inputs.g;
end


tau_y_lower_limit = sim_results.inputs.tau_y_lower_limit;
tau_y_upper_limit = sim_results.inputs.tau_y_upper_limit;

% Define the albedo of the lower boundary in our model
albedo_maxTau = sim_results.inputs.albedo_maxTau;




binEdges = sim_results.F_norm.binEdges;



% --------------------------------------------------------
% ------- PLOT TWO STREAM THEORY F_UP AND F_DOWN ---------
% --------------------------------------------------------


% K is defined in Bohren and Clothiaux (eq. 5.70)
K = sqrt((1 - ssa)*(1 - g*ssa));
% Define the reflectivity at the top of our layer, the photons that
% scatter out the cloud top
R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

% Define the constants
A = (R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


% lets plot some range of tau
tau = linspace(tau_y_lower_limit,tau_y_upper_limit,100);

C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];


s2 = subplot(3,1,2);

plot(photon_fraction_up(tau),tau,'Color',C1,'LineStyle','-','LineWidth',LW)
hold on;
plot(photon_fraction_down(tau),tau, 'Color',C2, 'LineStyle','-','LineWidth',LW)


% -------------------------------------------------------------
% ------- PLOT MONTE CARLO SIMULATION F_UP AND F_DOWN ---------
% -------------------------------------------------------------

plot(sim_results.F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle',':','LineWidth',LW)
plot(sim_results.F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle',':','LineWidth',LW)
grid on; grid minor
set(gca,'YDir','reverse')
ylabel('Optical Depth','Interpreter','latex', 'FontSize',17);

ylim([tau_y_lower_limit, tau_y_upper_limit])
xlim([0, 1.2])

% remove x axis tick label
set(gca,'XTickLabel',[])


% Set the yaxis tick labels
set(gca,'YTick', (0:1.3:6.5))
set(gca,'YTickLabel',(0:1.3:6.5))

set(gca,'FontSize',12)

%%



% -------------------------------------------------
% ***** LOAD RESULTS FOR OPTICAL DEPTH OF 48 *****
% -------------------------------------------------

sim_results = load('./Monte_Carlo_Simulation_Results/2D_MC_02-Dec-2022_Wavelength_550_N-Photons_1000000_N-Layers_1_Tau0_48_SZA_0.mat');

% -------------------------
% ***** Unpack Inputs *****
% -------------------------

if sim_results.inputs.mie.integrate_over_size_distribution==true
    ssa = sim_results.inputs.ssa_avg;
    g = sim_results.inputs.g_avg;
else
    ssa = sim_results.inputs.ssa;
    g = sim_results.inputs.g;
end


tau_y_lower_limit = sim_results.inputs.tau_y_lower_limit;
tau_y_upper_limit = sim_results.inputs.tau_y_upper_limit;

% Define the albedo of the lower boundary in our model
albedo_maxTau = sim_results.inputs.albedo_maxTau;




binEdges = sim_results.F_norm.binEdges;



% --------------------------------------------------------
% ------- PLOT TWO STREAM THEORY F_UP AND F_DOWN ---------
% --------------------------------------------------------


% K is defined in Bohren and Clothiaux (eq. 5.70)
K = sqrt((1 - ssa)*(1 - g*ssa));
% Define the reflectivity at the top of our layer, the photons that
% scatter out the cloud top
R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

% Define the constants
A = (R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_y_upper_limit)/...
    (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


% lets plot some range of tau
tau = linspace(tau_y_lower_limit,tau_y_upper_limit,100);

C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];


s3 = subplot(3,1,3);


plot(photon_fraction_up(tau),tau,'Color',C1,'LineStyle','-','LineWidth',LW)
hold on;
plot(photon_fraction_down(tau),tau, 'Color',C2, 'LineStyle','-','LineWidth',LW)


% -------------------------------------------------------------
% ------- PLOT MONTE CARLO SIMULATION F_UP AND F_DOWN ---------
% -------------------------------------------------------------

plot(sim_results.F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle',':','LineWidth',LW)
plot(sim_results.F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle',':','LineWidth',LW)
grid on; grid minor
set(gca,'YDir','reverse')
ylabel('Optical Depth','Interpreter','latex', 'FontSize',17);
xlabel('Normalized Irradiance $F/F_0$','Interpreter','latex', 'FontSize', 17);

ylim([tau_y_lower_limit, tau_y_upper_limit])
xlim([0, 1.2])


% Set the yaxis tick labels
set(gca,'YTick', (0:8:48))
set(gca,'YTickLabel',(0:8:48))

set(gca,'FontSize',12)

%%  ARRANGE THE FINAL POSITION OF EACH PLOT


% Now shift the positions of the subfigures

s1.Position = s1.Position.*[1 1 1 1.25];

s2.Position = s2.Position.*[1 1 1 1.25];

s3.Position = s3.Position.*[1 1 1 1.25];



