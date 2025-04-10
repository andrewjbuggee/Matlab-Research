%% Recreate Bohren and Clothiaux figure 6.10 with 3D monte Carlo

% By Andrew John Buggee


%% run a 3D Monte Carlo simulation for an optical thickness of 48 with overhead sun, ssa = 1, g = 0.85


clear variables


% ----- Define the boundaries of the medium ------
% comment out lines below if you want the medium to extend forever along
% that direction

inputs.tau_z_lower_limit = 0;
inputs.tau_z_upper_limit = 48;



% --------------------------------------------------


% define the solar zenith angle
% This is the angle of the incident radiation with respect to the medium
% normal direction
inputs.solar_zenith_angle = 0;                  % deg from zenith

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;

% Define the number of photons to inject into the medium
inputs.N_photons = 1e7;


% ----- Do you want to create a non-linear droplet profile? -----
inputs.createDropletProfile = false;

inputs.re = 10;

% define the wavelength
inputs.wavelength = 500;          % nanometers

% do you want to compute average ssa and g at each cloud layer?
% if so, the code will create a distribution of droplet sizes at each layer
% with the re value defining the modal radius
inputs.mie.integrate_over_size_distribution = false;



% Define the number of layers and the boundaries values for each tau
% layer. Properties are homogeneous in the X and Y plane, but vary in
% the Z plane

% Define the number of layers within the medium that differ
inputs.N_layers = 1;

% This options creates a simple cloud with a linear droplet profile
% or a homogenous cloud with a single radii
inputs.layerRadii = linspace(inputs.re,inputs.re, inputs.N_layers);      % radius of spheres in each layer


% Define the layer boundaries given the number of layers and the boundaries
% of the entire medium
inputs.layerBoundaries = linspace(inputs.tau_z_lower_limit, inputs.tau_z_upper_limit, inputs.N_layers +1);





% --- OVERRIDE SCATTERING PARAMETERS ----
inputs.g = 0.85;
inputs.g_avg = inputs.g;

inputs.ssa = 1;
inputs.ssa_avg = inputs.ssa;
% ----------------------------------------



% ------- Without Live Plotting ---------
tic
[F_norm, final_state, photon_tracking, inputs] = threeD_monteCarlo(inputs);
toc

%% Plot the upward and downward flux within the cloud

% define the center of each bin
bin_center = F_norm.binEdges(1:end-1) + diff(F_norm.binEdges)/2;

figure;

plot(F_norm.down, bin_center, ':', 'markersize', 10, 'Color', mySavedColors(1, 'fixed'))
hold on
plot(F_norm.up, bin_center, ':', 'markersize', 10, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor

xlabel('Normalized Irradiance', 'Interpreter','latex')
ylabel('Optical Depth', 'Interpreter','latex')

legend('$F_{\downarrow}$', '$F_{\uparrow}$', 'interpreter', 'latex', 'location', 'best')

title(['$\varpi = $', num2str(inputs.ssa_avg),',   $g = $', num2str(inputs.g_avg),...
    ',   $N_{photons} = 10^{', num2str(log10(inputs.N_photons)), '}$'], 'Interpreter','latex')

set(gcf, 'Position', [0 0 1300 750])

%% Plot the two stream analytical solution

% --------------------------------------------------------
% ------- PLOT TWO STREAM THEORY F_UP AND F_DOWN ---------
% --------------------------------------------------------




% K is defined in Bohren and Clothiaux (eq. 5.70)
K = sqrt((1 - inputs.ssa_avg)*(1 - inputs.g_avg*inputs.ssa_avg));
% Define the reflectivity at the top of our layer, the photons that
% scatter out the cloud top
R_inf = (sqrt(1-inputs.ssa_avg*inputs.g_avg) - sqrt(1 - inputs.ssa_avg))/...
    (sqrt(1-inputs.ssa_avg*inputs.g_avg) + sqrt(1 - inputs.ssa_avg));

% Define the constants
A = (R_inf - inputs.albedo_maxTau)*exp(-K*inputs.tau_z_upper_limit)/...
    (R_inf*(R_inf - inputs.albedo_maxTau)*exp(-K*inputs.tau_z_upper_limit) -...
    (1 - inputs.albedo_maxTau*R_inf)*exp(K*inputs.tau_z_upper_limit));

B = -(1 - R_inf*inputs.albedo_maxTau)*exp(K*inputs.tau_z_upper_limit)/...
    (R_inf*(R_inf - inputs.albedo_maxTau)*exp(-K*inputs.tau_z_upper_limit) -...
    (1 - inputs.albedo_maxTau*R_inf)*exp(K*inputs.tau_z_upper_limit));

photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


% plot!
LW = 3;         % linewidth
hold on
% plot(photon_fraction_up(bin_center),bin_center,'Color', mySavedColors(2, 'fixed'),...
%     'LineStyle','-','LineWidth',LW)
% hold on;
% plot(photon_fraction_down(bin_center),bin_center, 'Color', mySavedColors(1, 'fixed'),...
%     'LineStyle','-','LineWidth',LW)


[R, T, ~] = two_stream_RT(bin_center, linspace(inputs.ssa, inputs.ssa, length(bin_center)), ...
    linspace(inputs.g, inputs.g, length(bin_center)), linspace(inputs.albedo_maxTau,...
    inputs.albedo_maxTau, length(bin_center)));

plot(R,bin_center,'Color', mySavedColors(2, 'fixed'),...
    'LineStyle','-','LineWidth',LW)
hold on;
plot(T,bin_center, 'Color', mySavedColors(1, 'fixed'),...
    'LineStyle','-','LineWidth',LW)



legend('$F_{\downarrow}$ Monte Carlo', '$F_{\uparrow}$ Monte Carlo',...
    '$F_{\downarrow}$ 2 stream theory', '$F_{\uparrow}$ 2 stream theory',...
    'interpreter', 'latex', 'location', 'best')

