%% Recreate Bohren and Clothiaux figure 6.8 with 3D monte Carlo

% By Andrew John Buggee


%% First let's compute the 3D Monte Carlo transmittance for optical depths ranging from 1 to 50


clear variables


% ----- Define the boundaries of the medium ------
% comment out lines below if you want the medium to extend forever along
% that direction

tau_z_total = [1:10, 15:5:50];

inputs.tau_z_lower_limit = 0;



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








% --- OVERRIDE SCATTERING PARAMETERS ----
inputs.g = 0.85;
inputs.g_avg = 0.85;

inputs.ssa = 1;
inputs.ssa_avg = 1;
% ----------------------------------------


transmissivity_sza0 = zeros(1, length(tau_z_total));


for tt = 1:length(tau_z_total)

    disp([newline, 'tt = ', num2str(tt), '...', newline])

    % define the z boundary
    inputs.tau_z_upper_limit = tau_z_total(tt);

    % Define the layer boundaries given the number of layers and the boundaries
    % of the entire medium
    inputs.layerBoundaries = linspace(inputs.tau_z_lower_limit, inputs.tau_z_upper_limit, inputs.N_layers +1);

    % ------- Without Live Plotting ---------
    [~, final_state, ~, inputs] = threeD_monteCarlo(inputs);

    % store the transmissivity
    transmissivity_sza0(tt) = final_state.transmittance;

end

%%

clear inputs final_state, 


% ----- Define the boundaries of the medium ------
% comment out lines below if you want the medium to extend forever along
% that direction

tau_z_total = [1:10, 15:5:30];

inputs.tau_z_lower_limit = 0;



% --------------------------------------------------


% define the solar zenith angle
% This is the angle of the incident radiation with respect to the medium
% normal direction
inputs.solar_zenith_angle = 60;                  % deg from zenith

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;

% Define the number of photons to inject into the medium
inputs.N_photons = 1e6;


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








% --- OVERRIDE SCATTERING PARAMETERS ----
inputs.g = 0.85;
inputs.g_avg = 0.85;

inputs.ssa = 1;
inputs.ssa_avg = 1;
% ----------------------------------------


transmissivity_sza60 = zeros(1, length(tau_z_total));


for tt = 1:length(tau_z_total)

    disp([newline, 'tt = ', num2str(tt), '...', newline])

    % define the z boundary
    inputs.tau_z_upper_limit = tau_z_total(tt);

    % Define the layer boundaries given the number of layers and the boundaries
    % of the entire medium
    inputs.layerBoundaries = linspace(inputs.tau_z_lower_limit, inputs.tau_z_upper_limit, inputs.N_layers +1);

    % ------- Without Live Plotting ---------
    [~, final_state, ~, inputs] = threeD_monteCarlo(inputs);

    % store the transmissivity
    transmissivity_sza60(tt) = final_state.transmittance;

end



%% include the results with solar zenith angle of 60 degrees

hold on 

plot(tau_z_total, transmissivity_sza0, '.-', 'markersize', 25, 'LineWidth',1)

legend('$SZA = 0\deg$', '$SZA = 60\deg$', 'interpreter', 'latex', 'location', 'best')

title('$\varpi = 1$,  $g = 0.85$', 'Interpreter','latex')

%% Plot the two stream analytical solution

R = zeros(1, length(tau_z_total));
T = zeros(1, length(tau_z_total));
A = zeros(1, length(tau_z_total));

for tt = 1:length(tau_z_total)

    [R(tt), T(tt), A(tt)] = two_stream_RT(tau_z_total(tt), inputs.ssa_avg, inputs.g_avg, inputs.albedo_maxTau);

end

% add transmissivity to the plot
hold on 

plot(tau_z_total, T, '.-', 'markersize', 25, 'LineWidth',1)

legend('$SZA = 0\deg$', '$SZA = 60\deg$', '2 Stream Theory', 'interpreter', 'latex', 'location', 'best')

%% Plot the exponential attenuation using mie calcs


% compute exponential attenuation
tau_exp = 0:0.1:3;
T_exp = exp(-tau_exp);

% plot exponential attenuation on top of the other results
hold on 

plot(tau_exp, T_exp, '.-', 'markersize', 25, 'LineWidth',1)

legend('$SZA = 0\deg$', '$SZA = 60\deg$', '2 Stream Theory', 'Exponential Attenuation',...
    'interpreter', 'latex', 'location', 'best')

