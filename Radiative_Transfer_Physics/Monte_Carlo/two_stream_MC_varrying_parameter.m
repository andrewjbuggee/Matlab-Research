%% loop through 2-stream monte carlo with varrying parameters

% By Andrew John Buggee

%% Define inputs

clear variables


% Define the boundaries of the medium
inputs.tau_lower_limit = 0;
total_tau = logspace(-1, 2, 20);

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;

% Layers can differ by radius, by material, or both. If the layers differ
% by radius, create a vector describing each layer radius from top to
% bottom (tau = 0 to tau = tau0). If the layers differ by index of
% refraction, create a vector describing each layer indec of refraction
% from top to bottom. If both are true, create a cell array for both the
% changing radii and the changing index of refraction

inputs.layerRadii = 10;      % microns - radius of the spheres in each layer

% Define the number of layers within the medium that differ
inputs.N_layers = length(inputs.layerRadii);


% Define the number of photons to inject into the medium
inputs.N_photons = 1e6;

% Do you want to compute the internal upwelling and downwelling fluxes?
inputs.compute_internal_fluxes = false;


%%  MIE CALCULATIONS

% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

% define the wavelength
% The wavelength input is defined as follows:
% [wavelength_start, wavelength_end, wavelength_step].
inputs.mie.wavelength = [1600, 1600, 0];          % nanometers

% The first entry belows describes the type of droplet distribution
% that should be used. The second describes the distribution width. If
% running a mono-dispersed calculation, no entry for distribution width is
% required.
inputs.mie.distribution = {'mono', []};           % droplet distribution

% What mie code should we use to compute the scattering properties?
inputs.mie.mie_program = 'MIEV0';               % type of mie algorithm to run

% Do you want a long or short error file?
inputs.mie.err_msg_str = 'verbose';

% What is the index of refraction of the scatterer? If there is more than
% one scatterer, enter multiple values for the index of refraction
inputs.mie.indexOfRefraction = 'water';
%inputs.mie.indexOfRefraction = 1.33 + 9.32e-5i;

% Define the size of the scatterer and its scattering properties
% Assuming a pure homogenous medium composed of a single substance.
% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
if inputs.N_layers==1
    inputs.mie.radius = [inputs.layerRadii, inputs.layerRadii, 0];    % microns
else
    % define the min, max, and step. Record the vector because these are
    % the exact values used in the mie calculations

    % min value
    inputs.mie.radius(1) = min(inputs.layerRadii);
    % max value
    inputs.mie.radius(2) = max(inputs.layerRadii);
    % step value
    inputs.mie.radius(3) = abs(inputs.layerRadii(1) - inputs.layerRadii(2));    % microns

end




% Create a mie file
[input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,...
    inputs.mie.radius,inputs.mie.wavelength,inputs.mie.distribution, inputs.mie.err_msg_str, 1);

% run the mie file
runMIE(mie_folder,input_filename,output_filename);

% Read the output of the mie file
[ds,~,~] = readMIE(mie_folder,output_filename);

% --------------------------------------------------
% Outputs vary by wavelength along the row dimension
% Outputs vary by radii along the column dimension
% --------------------------------------------------

% Define the single scattering albedo
inputs.ssa = ds.ssa;

% Define the asymmetry parameter
inputs.g = ds.asymParam;





%% Run 2 stream 1D monte carlo code


% --------------------------------------
% ----- Override ssa and g values ------
% --------------------------------------
ssa = [1];
inputs.g = 0.85;
% --------------------------------------

R = zeros(length(ssa), length(total_tau));
T = zeros(length(ssa), length(total_tau));
A = zeros(length(ssa), length(total_tau));


legend_str = cell(1, length(ssa));

tic
for ss = 1:length(ssa)



    inputs.ssa = ssa(ss);

    legend_str{ss} = ['$\varpi = $', num2str(inputs.ssa)];

    for nn = 1:length(total_tau)

        disp(['[ss,nn] = [' num2str(ss),',',num2str(nn),'] ...'])

        inputs.tau_upper_limit = total_tau(nn);

        % Define the layer boundaries given the number of layers and the boundaries
        % of the entire medium
        inputs.layerBoundaries = linspace(inputs.tau_lower_limit, inputs.tau_upper_limit, inputs.N_layers +1);


        % ------- Without Live Plotting ---------

        [F_norm, final_state, photon_tracking, inputs] = twoStream_monteCarlo(inputs);


        R(ss, nn) = final_state.scatter_out_top/inputs.N_photons;
        T(ss, nn) = final_state.scatter_out_bottom/inputs.N_photons;
        A(ss, nn) = final_state.absorbed/inputs.N_photons;






    end
end
toc

%% Plot

figure;

plotx(total_tau, R, '.-')

% Label axes
xlabel('Optical Thickness');
ylabel('Reflectivity');

% Improve plot aesthetics
grid on; grid minor
title('1D Monte Carlo with Absorption', ['g=', num2str(inputs.g)]);
legend(legend_str, 'Location','best', 'Interpreter','latex')

% plot the theoretical values on top 
R_theory = zeros(length(ssa), length(total_tau));

for ss = 1:length(ssa)

    R_theory(ss,:) = two_stream_RT(total_tau, linspace(ssa(ss), ssa(ss), length(total_tau)),...
        linspace(inputs.g, inputs.g, length(total_tau)), 0);

end

hold on
semilogx(total_tau, R_theory, '--', 'Color', 'k')

%% Make a plot of just the two-stream theory

% vary optical thickness and single scattering albedo



% plot the theoretical values on top 
R_theory = zeros(length(ssa), length(total_tau));

for ss = 1:length(ssa)

    R_theory(ss,:) = two_stream_RT(total_tau, linspace(ssa(ss), ssa(ss), length(total_tau)),...
        linspace(inputs.g, inputs.g, length(total_tau)), 0);

end

figure;

semilogx(total_tau, R_theory, '-')

% Label axes
xlabel('Optical Thickness');
ylabel('Reflectivity');

% Improve plot aesthetics
grid on; grid minor
title('1D Monte Carlo with Absorption', ['g=', num2str(inputs.g)]);
legend(legend_str, 'Location','best', 'Interpreter','latex')

%% The reflectivity of a pile of parallel plates

% Bohren and Clothiaux Chapter 5.1

R = total_tau./(1 + total_tau);


