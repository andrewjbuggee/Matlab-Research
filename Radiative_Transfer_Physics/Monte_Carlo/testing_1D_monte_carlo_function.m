%% Testing 2-stream monte carlo

% By Andrew John Buggee

%% Define inputs

clear variables


% Define the boundaries of the medium
inputs.tau_lower_limit = 0;
inputs.tau_upper_limit = 15;

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

% Define the layer boundaries given the number of layers and the boundaries
% of the entire medium
inputs.layerBoundaries = linspace(inputs.tau_lower_limit, inputs.tau_upper_limit, inputs.N_layers +1);

% Define the number of photons to inject into the medium
inputs.N_photons = 1e5;


%%  MIE CALCULATIONS

% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

% define the wavelength
% The wavelength input is defined as follows:
% [wavelength_start, wavelength_end, wavelength_step].
inputs.mie.wavelength = [550, 550, 0];          % nanometers

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

% --------------------------------------
% ----- Override ssa and g values ------
% --------------------------------------
inputs.ssa = 1;
inputs.g = 0.5;
% --------------------------------------



%% Run 2 stream 1D monte carlo code

% ------- Without Live Plotting ---------
[F_norm, final_state, photon_tracking, inputs] = twoStream_monteCarlo(inputs);


% ------- With Live Plotting ---------

% do you want to save the live figure for as a video file?
% inputs.saveLivePlot_asVideo = false;
% [F_norm, final_state, photon_tracking, inputs] = twoStream_monteCarlo_with_livePlot(inputs);

%% Compare the monte carlo solution with the analytical solution

plot_2strm_monteCarlo_with_analytical(inputs,F_norm);


%% Lets plot the max depth reached by each photon normalized by the total number of photons

plot_probability_maxDepth(inputs, photon_tracking, 'probability')


%% Let's plot the conditional probability of an absorbed photon reach a max depth of tau

%plot_probability_absorbedDepth(inputs, final_state, photon_tracking, 'probability')

%% Let's plot the conditional probability of a photon that scattered out the cloud top
% reaching a max depth of tau

%plot_probability_scatterOutTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Let's plot two conditional probabilities on the same plot
% Plot the probability of a photon reaching a max depth of tau if it was
% scattered out the cloud top, and plot the probability of a photon being
% absorbed at a depth tau given that it was absorbed.

plot_probability_absANDscatTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Plot a bar chart showing the probability of each final state

plot_probability_finalStates(final_state,inputs)

