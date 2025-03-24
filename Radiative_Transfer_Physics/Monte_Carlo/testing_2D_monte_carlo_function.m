%% Testing 2D monte carlo

% By Andrew John Buggee

%% ------ Define core inputs -----

clear variables


% Define the boundaries of the medium
inputs.tau_y_lower_limit = 0;
inputs.tau_y_upper_limit = 8; 

% define the solar zenith angle
% This is the angle of the incident radiation with respect to the medium
% normal direction
inputs.solar_zenith_angle = 49.45;                  % deg from zenith
%inputs.solar_zenith_angle = 27;                  % deg from zenith

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;

% Define the number of photons to inject into the medium
inputs.N_photons = 1e7;


% ----- Do you want to create a non-linear droplet profile? -----
inputs.createDropletProfile = true;

% --- if true.... ---
% Physical constraint that shapes the droplet profile
%inputs.dropletProfile.constraint = 'adiabatic';
inputs.dropletProfile.constraint = 'linear_with_z';
% Define the radius value at cloud top and cloud bottom
inputs.dropletProfile.r_top = 12;            % microns
inputs.dropletProfile.r_bottom = 5;          % microns

% --- else ---
inputs.re = 10;



% define the wavelength
inputs.wavelength = 1600;          % nanometers

% do you want to compute average ssa and g at each cloud layer?
% if so, the code will create a distribution of droplet sizes at each layer
% with the re value defining the modal radius
inputs.mie.integrate_over_size_distribution = true;

% --- if true... ---
% Define the type of size distribution
inputs.mie.size_dist = 'gamma';

% Define the distribution variance, depending on the distribution type used
% Has to be the same length as the numer of layers in our medium
inputs.size_distribution_var = 7;           % Typically value for liquid water clouds





%%
% ----------------------------------------------------------------------
% ------------- DEFINE EACH LAYER'S DROPLET RADIUS ---------------------
% ----------------------------------------------------------------------

% Layers can differ by radius, by material, or both. If the layers differ
% by radius, create a vector describing each layer radius from top to
% bottom (tau = 0 to tau = tau0). If the layers differ by index of
% refraction, create a vector describing each layer indec of refraction
% from top to bottom. If both are true, create a cell array for both the
% changing radii and the changing index of refraction






if inputs.createDropletProfile==false
    
    % Define the number of layers and the boundaries values for each tau
    % layer

    % Define the number of layers within the medium that differ
    inputs.N_layers = 1;

    % This options creates a simple cloud with a linear droplet profile
    % or a homogenous cloud with a single radii
    inputs.layerRadii = linspace(inputs.re,inputs.re, inputs.N_layers);      % radius of spheres in each layer


    

    % Define the layer boundaries given the number of layers and the boundaries
    % of the entire medium
    inputs.layerBoundaries = linspace(inputs.tau_y_lower_limit, inputs.tau_y_upper_limit, inputs.N_layers +1);


else

    % --------------------------------------------
    % --------- create droplet profile -----------
    % --------------------------------------------


    % define the number of layers to model within the cloud
    inputs.dropletProfile.N_layers = 100;


    % Define the boundaries of each tau layer
    % Define the layer boundaries given the number of layers and the boundaries
    % of the entire medium
    inputs.dropletProfile.layerBoundaries = linspace(inputs.tau_y_lower_limit, inputs.tau_y_upper_limit, inputs.dropletProfile.N_layers +1);

    % Define the optical depth vector that defines the mid point of each
    % layer
    inputs.dropletProfile.tau_layer_mid_points = inputs.dropletProfile.layerBoundaries(1:end-1) + diff(inputs.dropletProfile.layerBoundaries)/2;

    % tell the code if the vertical dimension is defined as altitude or
    % optical depth
    inputs.dropletProfile.independent_variable = 'optical_depth';                    % string that tells the code which independent variable we used

    % Compute the droplet profile
    inputs.dropletProfile.re = create_droplet_profile2([inputs.dropletProfile.r_top, inputs.dropletProfile.r_bottom],...
        inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.independent_variable,...
        inputs.dropletProfile.constraint);     % microns - effective radius vector

    % Now set the layerRadii with the computed droplet profile
    inputs.layerRadii = inputs.dropletProfile.re;


    % Define the number of layers and the layer boundaries at the base
    % level of the structure

    % Define the number of layers within the medium that differ
    inputs.N_layers = length(inputs.layerRadii);

    % Define the layer boundaries given the number of layers and the boundaries
    % of the entire medium
    inputs.layerBoundaries = inputs.dropletProfile.layerBoundaries;

end


% ----------------------------------------------------------------------
% There HAS to be even spacing for the mie calcualtion to take place!




%% MIE CALCULATIONS

% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

% The wavelength input is defined as follows:
% [wavelength_start, wavelength_end, wavelength_step]
inputs.mie.wavelength = [inputs.wavelength, inputs.wavelength, 0];          % nanometers

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



    % Create a mie file
    [input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,...
        inputs.mie.radius, inputs.mie.wavelength, inputs.mie.distribution, inputs.mie.err_msg_str, 1);

    % run the mie file
    [~] = runMIE(mie_folder,input_filename,output_filename);


    % Read the output of the mie file
    % --------------------------------------------------
    % Outputs vary by wavelength along the row dimension
    % Outputs vary by radii along the column dimension
    % --------------------------------------------------
    [ds,~,~] = readMIE(mie_folder,output_filename);


    % Define the asymmetry parameter
    inputs.g = ds.asymParam;

    % define the single scattering albedo
    inputs.ssa = ds.ssa;

    % define the extinction efficiency
    inputs.Qe = ds.Qext;



elseif inputs.N_layers>1 && inputs.createDropletProfile==false

    % If both of these conditions are met, we have a simple droplet profile
    % that increases linearly with tau

    % define the min, max, and step. Record the vector because these are
    % the exact values used in the mie calculations

    % Define the radius vector for our mie calculation
    inputs.mie.radiusVector = sort(inputs.layerRadii);
    % min value
    inputs.mie.radius(1) = min(inputs.layerRadii);
    % max value
    inputs.mie.radius(2) = max(inputs.layerRadii);
    % step value
    inputs.mie.radius(3) = inputs.mie.radiusVector(2) - inputs.mie.radiusVector(1);    % microns




    % Create a mie file
    [input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,...
        inputs.mie.radius,inputs.mie.wavelength,inputs.mie.distribution, inputs.mie.err_msg_str, 1);

    % run the mie file
    [~] = runMIE(mie_folder,input_filename,output_filename);


    % Read the output of the mie file
    % --------------------------------------------------
    % Outputs vary by wavelength along the row dimension
    % Outputs vary by radii along the column dimension
    % --------------------------------------------------
    [ds,~,~] = readMIE(mie_folder,output_filename);


    % ------------------------------------------------------------
    % Organize the scattering properties from medium top to bottom
    % Use the definition of the each layer's radii
    % This vector is define from layer top to layer bottom
    % ------------------------------------------------------------

    inputs.g = zeros(1, length(inputs.layerRadii));
    inputs.ssa = zeros(1, length(inputs.layerRadii));
    inputs.Qe = zeros(1, length(inputs.layerRadii));
    for LL = 1:inputs.N_layers

        % define the asymmetry parameter
        [~,index] = min(abs((inputs.mie.radius(1):inputs.mie.radius(3):inputs.mie.radius(2))-inputs.layerRadii(LL)));
        inputs.g(:,LL) = ds.asymParam(:,index);

        % define the single scattering albedo
        inputs.ssa(:,LL) = ds.ssa(:,index);

        % define the extinction efficiency
        inputs.Qe(:,LL) = ds.Qext(:,index);

    end



elseif inputs.N_layers>1 && inputs.createDropletProfile==true

    % If both of these conditions are met, than we have a droplet profile
    % that is non-linear. Therefore we have to run a mie calcualtion for
    % each radius value and collect the values ourselves in a loop

    inputs.g = zeros(1, length(inputs.layerRadii));
    inputs.ssa = zeros(1, length(inputs.layerRadii));
    inputs.Qe = zeros(1, length(inputs.layerRadii));

    for nn = 1:inputs.N_layers

        % Define the radius value
        inputs.mie.radius(1) = inputs.layerRadii(nn);           % microns


        % Create a mie file
        [input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,...
            inputs.mie.radius,inputs.mie.wavelength,inputs.mie.distribution, inputs.mie.err_msg_str, nn);

        % run the mie file
        [~] = runMIE(mie_folder,input_filename,output_filename);


        % Read the output of the mie file
        % --------------------------------------------------
        % Outputs vary by wavelength along the row dimension
        % Outputs vary by radii along the column dimension
        % --------------------------------------------------
        [ds,~,~] = readMIE(mie_folder,output_filename);


        % define the asymmetry parameter
        inputs.g(nn) = ds.asymParam;

        % define the single scattering albedo
        inputs.ssa(nn) = ds.ssa;

        % define the extinction efficiency
        inputs.Qe(nn) = ds.Qext;



    end


end



% We don't need the rest of the mie computations, so lets delete them
clear ds



% Do you want to integrate over a size distribution?

if inputs.mie.integrate_over_size_distribution==true

    inputs.mie.dist_var = linspace(inputs.size_distribution_var, inputs.size_distribution_var, inputs.N_layers);           % Typically value for liquid water clouds

    % Compute the average value for the single scattering albedo over a size
    % distribution
    [inputs.ssa_avg, inputs.Qe_avg, inputs.g_avg] = average_mie_over_size_distribution(inputs.layerRadii, inputs.mie.dist_var,...
        inputs.mie.wavelength(1),inputs.mie.indexOfRefraction, inputs.mie.size_dist, 1);


end


%% Run 2 stream 2D monte carlo code

tic

% ------- Without Live Plotting ---------
[F_norm, final_state, photon_tracking, inputs] = twoD_monteCarlo(inputs);

% ---------- With Live Plotting ---------
% [F_norm, final_state, photon_tracking, inputs] = twoD_monteCarlo_withLivePlot(inputs);

toc


%% Compare the 2D monte carlo solution with the 2-stream analytical solution

plot_2strm_2D_monteCarlo(inputs,F_norm);


%% Do you want to save you results?


if strcmp(whatComputer, 'anbu8374')
    % save in the following folder
    inputs.folder_name_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results/'];
    
elseif strcmp(whatComputer, 'andrewbuggee')==true

    % save in the following folder
    inputs.folder_name_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results/'];

end

if inputs.N_layers==1

    save([inputs.folder_name_2save, '2D_MC_',char(datetime('today')),'_Wavelength_',num2str(inputs.mie.wavelength(1)),...
        '_N-Photons_',num2str(inputs.N_photons),'_N-Layers_',num2str(inputs.N_layers),...
        '_Tau0_',num2str(inputs.tau_y_upper_limit),'_r_e_',num2str(inputs.dropletProfile.re),...
        '_SZA_',num2str(inputs.solar_zenith_angle),'.mat'],...
        "inputs","F_norm", "final_state", "photon_tracking");

elseif inputs.N_layers>1 && inputs.createDropletProfile==true

    save([inputs.folder_name_2save, '2D_MC_',char(datetime('today')),'_Wavelength_',num2str(inputs.mie.wavelength(1)),...
        '_N-Photons_',num2str(inputs.N_photons),'_N-Layers_',num2str(inputs.N_layers),...
        '_Tau0_',num2str(inputs.tau_y_upper_limit),'_r_top_',num2str(inputs.dropletProfile.r_top),...
        '_r_bot_',num2str(inputs.dropletProfile.r_bottom),'_SZA_',num2str(inputs.solar_zenith_angle),'.mat'],...
        "inputs","F_norm", "final_state", "photon_tracking");
end





%% Find the pdf of the number of scattering events

figure; 
histogram(photon_tracking.number_of_scattering_events, 'Normalization', 'pdf') 
grid on; grid minor;
ylabel('PDF', 'Interpreter', 'latex')
xlabel('Number of Scattering Events', 'Interpreter', 'latex')
title(['$\tau_c = $', num2str(inputs.tau_y_upper_limit)], 'Interpreter', 'latex')
set(gcf, 'Position', [0 0 1000, 750])

% This is a discrete PDF! Values can only be integers
texBox_str = {['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\varpi$ = ',num2str(inputs.ssa)],...
    ['mean = ', num2str(mean(photon_tracking.number_of_scattering_events))],...
    ['median = ', num2str(median(photon_tracking.number_of_scattering_events))],...
    ['mode = ', num2str(mode(photon_tracking.number_of_scattering_events))]};

dim = [0.6 0.85 0 0];
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';

%% Lets plot the max depth reached by each photon normalized by the total number of photons

%plot_probability_maxDepth(inputs, photon_tracking, 'probability')


%% Let's plot the conditional probability of an absorbed photon reach a max depth of tau

%plot_probability_absorbedDepth(inputs, final_state, photon_tracking, 'probability')

%% Let's plot the conditional probability of a photon that scattered out the cloud top
% reaching a max depth of tau

plot_probability_scatterOutTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Let's plot two conditional probabilities on the same plot
% Plot the probability of a photon reaching a max depth of tau if it was
% scattered out the cloud top, and plot the probability of a photon being
% absorbed at a depth tau given that it was absorbed.

plot_probability_absANDscatTop_maxDepth(inputs, final_state, photon_tracking, 'pdf')


%% Plot a bar chart showing the probability of each final state

plot_probability_finalStates(final_state,inputs)

