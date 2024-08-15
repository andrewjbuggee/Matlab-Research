% --- Compute change in reflectance due to a change in the radius at cloud bottom for EMIT channels-----



% By Andrew J. Buggee
%%

function [measurement_change, change_in_r_bottom] = compute_reflectanceChange_due_to_rBottom_change(emit, state_vector,...
    measurement_estimate, inputs, pixels2use, pp)

% --- define the reflectance uncertinaty ---
reflectance_uncertainty = emit.reflectance.uncertainty(inputs.bands2run);


% --- compute the measurement change at a specific state vector ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% ---------------------------------------------------------
% ---- define the incremental change to each variable -----

change_in_r_bottom = 0.25:0.25:3;       % microns
 
% ---------------------------------------------------------




% ----------------------------------------------------------

% Set up a few constants for the water cloud
H = inputs.RT.cloudDepth;                                % km - geometric thickness of cloud
n_layers = inputs.RT.cloud_layers;                          % number of layers to model within cloud

% Cloud top
z_top = inputs.RT.cloudTop_height(pp);        % km -  cloud top height

z = linspace(z_top-H, z_top,n_layers);        % km - altitude above ground vector

indVar = 'altitude';                    % string that tells the code which independent variable we used

profile_type = inputs.model.profile.type; % type of water droplet profile
dist_str = inputs.RT.drop_distribution_str;                         % droplet distribution

% -- For now, lets assume this is constant --
dist_var = linspace(inputs.RT.drop_distribution_var, inputs.RT.drop_distribution_var,...
    inputs.RT.cloud_layers);              % distribution variance
vert_homogeneous_str = inputs.RT.vert_homogeneous_str;          % This tells the function whether of not to create a multi-layered cloud
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector.

% Tell the code to use a pre-computed mie table for the extinction
% efficiency, or to use the value of the extinction paradox -> Qe = 2
parameterization_str = inputs.RT.parameterization_str;

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = emit.radiance.wavelength(inputs.bands2run(1));    % nm - Wavelength used for cloud optical depth calculation

% Lets step through each model variable and compute the derivative
measurement_change = zeros(length(measurement_estimate), length(change_in_r_bottom));

% ----- Determine when we have a change in our measurement above the uncertainty -----


for xx = 1:length(change_in_r_bottom)
    
    tic 
    % add to the radius at cloud bottom
    new_r_bottom = r_bottom + change_in_r_bottom(xx);       % microns

    % ------------------------------------------------
    % create water cloud file with new droplet profile
    % ------------------------------------------------
    
    new_re = create_droplet_profile2([r_top, new_r_bottom], z, indVar, profile_type);     % microns - effective radius vector
    
    
    loop_var = 0;

    wc_filename = write_wc_file(new_re, tau_c, z_topBottom, wavelength_tau_c(1,1), dist_str,...
        dist_var, vert_homogeneous_str, parameterization_str, loop_var);
    
    
    % ----- Write an INP file --------
    names.inp = write_INP_file_4EMIT_Gauss_Newton(inputs, pixels2use, emit, wc_filename);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
    
    % ---- Run uvspec for the files created -----
    [new_measurement_estimate,~] = runReflectanceFunction_4EMIT_gaussNewton(names, inputs, emit.spec_response.value);
    
    measurement_change(:,xx) = new_measurement_estimate' - measurement_estimate;
    
    % --- check is all reflectance changes are above measurement uncertainty 
    if all(measurement_change(:,xx)>reflectance_uncertainty)==true
    
        % clear the rest of the zerod matrix
        measurement_change(:,(xx+1:length(change_in_r_bottom))) = [];
        break
    
    end

    disp([newline, 'Iteration: r_bot = ',num2str(xx),'/', num2str(length(change_in_r_bottom)), ...
        ', Time to run was: ', num2str(toc), ' sec', newline])
    
    
end




end