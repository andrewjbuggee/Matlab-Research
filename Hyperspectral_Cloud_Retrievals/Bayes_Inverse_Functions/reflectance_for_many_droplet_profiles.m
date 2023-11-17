%% Compute the Reflectance for a cloud with a droplet profile


% By Andrew John Buggee


%%

r_top = 5:15;       % microns
r_bot = 5:10;       % microns
tau_c = 5:2:35;     % optical depth


% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------

H = 0.5;                                % km - geometric thickness of cloud
n_layers = 5;                           % number of layers to model within cloud

z0 = 1;                                 % km - base height of cloud
z = linspace(z0, z0+H,n_layers);        % km - altitude above ground vector
indVar = 'altitude';                    % string that tells the code which independent variable we used
constraint = 'adiabatic';              % string that tells the code which physical constraint to use

wavelength_tau_c = 645;                 % nm
dist = 'mono';                         % droplet distribution
homogenous_str = 'non-homogeneous';     % This tells the function to create a multi-layered cloud
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector. 


% First lets fix the r_bot to be 5 microns and tau to be 10
for ii = 1:length(r_top)
    re(:,ii) = create_droplet_profile2([r_top, r_bottom(1)], z, indVar, constraint);     % microns - effective radius vector
    wc_filename(ii) = write_wc_file(re(ii), 10, z_topBottom, wavelength_tau_c(1,1), dist, homogenous_str);

end

% Now lets calculate the reflectance
% ----- Write an INP file --------
GN_names.inp = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename);
    
% now lets write the output names
    
GN_names.out = writeOutputNames(GN_names.inp);

% ---- Run uvspec for the files created -----
[measurement_estimate,~] = runReflectanceFunction_4gaussNewton(GN_names,INP_folderName,saveCalculations_fileName);



%% Now lets fix r_top and vary r_bottom





