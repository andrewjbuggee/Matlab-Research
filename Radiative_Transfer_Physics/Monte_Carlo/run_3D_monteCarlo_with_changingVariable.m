%% Run 3D monte carlo in a loop with a changing variable

% By Andrew John Buggee

%%

%% Set up inputs for 3D monte carlo, run program, plot results

% By Andrew John Buggee

%% ------ Define core inputs -----

clear variables


% ----- Define the boundaries of the medium ------
% comment out lines below if you want the medium to extend forever along
% that direction

% inputs.tau_y_lower_limit = 0;
% inputs.tau_y_upper_limit = 8;
%
% inputs.tau_x_lower_limit = 0;
% inputs.tau_x_upper_limit = 8;

inputs.tau_z_lower_limit = 0;
% inputs.tau_z_upper_limit = logspace(-1,2,25);
total_tau = logspace(-1,2,30);

% --------------------------------------------------


% define the solar zenith angle
% This is the angle of the incident radiation with respect to the medium
% normal direction
inputs.solar_zenith_angle = 0;                  % deg from zenith
%inputs.solar_zenith_angle = 27;                  % deg from zenith

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;

% Define the number of photons to inject into the medium
inputs.N_photons = 1e5;


% ----- Do you want to create a non-linear droplet profile? -----
inputs.createDropletProfile = false;

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
inputs.wavelength = 500;          % nanometers

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



% Do you want to compute the internal fluxes within the medium? If not, and
% you only care about total absorption, transmission and reflectance, set
% this flag to false
inputs.compute_internal_fluxes = false;





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


R = zeros(1, length(total_tau));
T = zeros(1, length(total_tau));
A = zeros(1, length(total_tau));


for nn = 1:length(total_tau)


    inputs.tau_z_upper_limit = total_tau(nn);



    if inputs.createDropletProfile==false

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


    %% Run 3D monte carlo code

    % --- OVERRIDE SCATTERING PARAMETERS ----
    inputs.g = 0.85;
    inputs.g_avg = inputs.g;

    inputs.ssa = 0.99;
    inputs.ssa_avg = inputs.ssa;
    % ----------------------------------------


    tic

    % ------- Without Live Plotting ---------
    [~, final_state, ~, inputs] = threeD_monteCarlo(inputs);

    toc

    % Keep just the reflectance, transmittance and absorptance
    R(nn) = final_state.reflectance;
    T(nn) = final_state.transmittance;
    A(nn) = final_state.absorptance;



end

%% Plot



figure;

plot(total_tau, R, '.-', 'Color', mySavedColors(1, 'fixed'), 'MarkerSize', 20, 'LineWidth', 1.5)
hold on

% Label axes
xlabel('Optical Thickness', "Interpreter","latex");


% Improve plot aesthetics
grid on; grid minor
title('3D Monte Carlo with Absorption ', ['$g=$', num2str(inputs.g_avg), ',   $\varpi = $',...
    num2str(inputs.ssa_avg)], 'Interpreter','latex', 'FontSize', 20);

% plot the transmissivity
plot(total_tau, T, '.-', 'Color', mySavedColors(2, 'fixed'), 'MarkerSize', 20, 'LineWidth', 1.5)

% plot the absoroptivity
plot(total_tau, A, '.-', 'Color', mySavedColors(3, 'fixed'), 'MarkerSize', 20, 'LineWidth', 1.5)

ylabel('$\%$', 'Interpreter','latex');

legend('$R_{3D}$', '$T_{3D}$', '$A_{3D}$',...
     'Location','best', 'Interpreter','latex')


% % ------ plot the 2-stream theoretical values ------
% 
% [R_theory, T_theory, A_theory] = two_stream_RT(total_tau, linspace(inputs.ssa_avg, inputs.ssa_avg, length(total_tau)),...
%     linspace(inputs.g_avg, inputs.g_avg, length(total_tau)), 0);
% 
% 
% hold on
% 
% semilogx(total_tau, R_theory, 'Color', mySavedColors(1, 'fixed'), 'LineWidth', 1);
% 
% semilogx(total_tau, T_theory, 'Color', mySavedColors(2, 'fixed'), 'LineWidth', 1);
% 
% semilogx(total_tau, A_theory, 'Color', mySavedColors(3, 'fixed'), 'LineWidth', 1)
% 
% 
% 
% % legend([legend_str, {'Analytical'}], 'Location','best', 'Interpreter','latex')
% legend('$R_{3D}$', '$T_{3D}$', '$A_{3D}$', '$R_{2-stream}$', '$T_{2-stream}$', '$A_{2-stream}$',...
%     'Location','best', 'Interpreter','latex')

