%% Write the INP files used by libRadtran to solve the radiative transfer equation

% By Andrew John Buggee
%%

function [] = write_INP_file(INP_folderpath, libRadtran_data_path, inputFileName, inputs,...
                        wavelengths, wc_filename, mc_basename, wc_modify_tau, waterVaporProfile_filename)



% ----------------- ******************** -------------------
% ------------------ Write the INP File --------------------
% ----------------- ******************** -------------------

% Open the old file for writing
fileID = fopen([INP_folderpath, inputFileName], 'w');

% Define which RTE solver to use
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
fprintf(fileID, formatSpec,'rte_solver', inputs.RT.rte_solver,' ', '# Radiative transfer equation solver');


if strcmp(inputs.RT.rte_solver, 'montecarlo')==true

    % Define the number of photons to use in the simulation
    % ------------------------------------------------
    formatSpec = '%s %u %5s %s \n';
    fprintf(fileID, formatSpec,'mc_photons', inputs.RT.mc.photons,' ', '# MYSTIC number of photons');


    if isfield(inputs.RT.mc, 'vroom')

        % Define whether or not to use VROOM
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'mc_vroom', inputs.RT.mc.vroom,' ', '# helps speed up calculations for particles with strong forward scattering');

    end

    if isfield(inputs.RT.mc, 'escape')

        % Do you wish to calculate the escape probabilities?
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'mc_escape', inputs.RT.mc.escape,' ', '# calculates radiances via escape probabilities - speeds up computation');

    end

    if isfield(inputs.RT.mc, 'backward')
        % Compute backwards ray-tracing
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'mc_backward', inputs.RT.mc.backward,' ', '# backward tracing of photons');

    end

    if exist('mc_basename','var')


        % Define the number of streams to keep track of when solving the equation
        % of radiative transfer
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'mc_basename', mc_basename,' ', '# file names and path for 3D mystic output');

    end


else


    % Define the number of streams to keep track of when solving the equation
    % of radiative transfer
    % ------------------------------------------------
    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'number_of_streams', inputs.RT.num_streams,' ', '# Number of streams');

end


% Use phase function correction?
% ------------------------------------------------
if isfield(inputs.RT, 'use_nakajima_phaseCorrection')
    if inputs.RT.use_nakajima_phaseCorrection==true
        % define the pahse correction to be true
        % ------------------------------------------------
        formatSpec = '%s %5s %s \n\n';
        fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
    end
end

% Define location of the data files
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
fprintf(fileID, formatSpec,'data_files_path', libRadtran_data_path, ' ', '# Location of libRadtran data files');


% Define the band model to use
% of radiative transfer
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
fprintf(fileID, formatSpec,'mol_abs_param', inputs.RT.band_parameterization,' ', '# Band parameterization');


% Define the location and filename of the atmopsheric profile to use
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
fprintf(fileID, formatSpec,'atmosphere_file ', [libRadtran_data_path, 'atmmod/', inputs.RT.atm_file],...
    ' ', '# Location of atmospheric profile');

% Define the location and filename of the extraterrestrial solar source
% ---------------------------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
fprintf(fileID, formatSpec,'source solar', inputs.RT.source_file, ' ', '# Bounds between 250 and 10000 nm');


% Define the day of the calendar year to account for Earth-Sun distance
% ---------------------------------------------------------------------
if isfield(inputs.RT, 'day_of_year')

    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'day_of_year', inputs.RT.day_of_year, ' ', '# accounts for changing Earth-Sun distance');

end



% Define the surface albedo
% ------------------------------------------------
formatSpec = '%s %f %5s %s \n\n';
fprintf(fileID, formatSpec,'albedo', inputs.RT.surface_albedo, ' ', '# Surface albedo of the ocean');


% Define the Water Cloud properties, if you want a cloud in your model
% --------------------------------------------------------------------
if inputs.RT.yesCloud==true

    % Define the water cloud file
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/', wc_filename], ' ', '# Location of water cloud file');
    %fprintf(fileID, formatSpec,'wc_file 1D', [libRadtran_data_path,'wc/', wc_filename{rr,tc}{1}], ' ', '# Location of water cloud file');



    % Define the technique or parameterization used to convert liquid cloud
    % properties of r_eff and LWC to optical depth
    % ----------------------------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'wc_properties', inputs.RT.wc_parameterization, ' ', '# optical properties parameterization technique');


    if isfield(inputs.RT, 'modify_wc_opticalDepth')==true && inputs.RT.modify_wc_opticalDepth==true

        % Manually set the water cloud optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'wc_modify tau set', wc_modify_tau, ' ', '# optical properties parameterization technique');

    end


end


if inputs.RT.monochromatic_calc==false

    % Define the wavelengths for which the equation of radiative transfer will
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'wavelength', wavelengths(1), wavelengths(2), ' ', '# Wavelength range');

else

    % Define the wavelengths for which the equation of radiative transfer will
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'wavelength', wavelengths, wavelengths, ' ', '# Wavelength range');

end

% Spline interpolate the calculated spectrum between wavelengths lambda_0
% and lambda_1 in steps of lambda_step, in nm.
% -------------------------------------------------------------------------
% formatSpec = '%s %f %f %5s %s \n\n';
% fprintf(fileID, formatSpec,'spline', inputs.RT.wavelength(ww, 1), inputs.RT.wavelength(ww, 2), ' ', '# Wavelength range');




if inputs.RT.use_coxMunk==true

    % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'brdf_cam u10', inputs.RT.wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

end


% ----- Do you want to specify specific cross-section models? -----
if inputs.RT.specify_cross_section_model==true

    if isfield(inputs.RT, 'crs_model_rayleigh')==true

        % Define which specific models you'd like to use
        % valid for rayleigh scattering, Ozone, NO2, and O4
        % -------------------------------------------------------------------------
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'crs_model rayleigh', inputs.RT.crs_model_rayleigh, ' ', '# Rayleigh scattering cross section');

    end

end





% Define the total column water vapor amount
% --------------------------------------------------------------------
if isfield(inputs.RT, 'modify_total_columnWaterVapor') && inputs.RT.modify_total_columnWaterVapor==true

    % If true, modify the amount of column water vapor
    % --------------------------------------------------------------
    formatSpec = '%s %f %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'mol_modify H2O ', inputs.RT.waterVapor_column, ' MM', ' ', '# Column water vapor amount');


end



% Alter the above cloud column water vapor amount
% --------------------------------------------------------------------
if isfield(inputs.RT, 'modify_aboveCloud_columnWaterVapor') && inputs.RT.modify_aboveCloud_columnWaterVapor==true

    % check to make sure the input for this setting exists
    if exist("waterVaporProfile_filename", "var")==false

        error([newline, 'No custom file defined for a water vapor density profile.', newline])

    end

    % If true, define the filename for the new custom column water vapor
    % density profile
    % ----------------------------------------------------------------
    formatSpec = '%s %s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'mol_file H2O ', waterVaporProfile_filename, ' cm_3', ' ', '# Custom water vapor profile');


end





% Define the concentration of carbon dioxide
% --------------------------------------------------------------------
if inputs.RT.modify_CO2==true

    % If true, modify the mixing ratio of carbon dioxide
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'mixing_ratio CO2 ', inputs.RT.CO2_mixing_ratio, ' ', '# ppm of CO2');


end


% Define the concentration of molecular nitrogen
% --------------------------------------------------------------------
if inputs.RT.modify_N2==true

    % If true, modify the mixing ratio of molecular nitrogen
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s %s \n\n';
    fprintf(fileID, formatSpec,'mol_modify N2 ', inputs.RT.N2_mixing_ratio, 'CM_2', ' ', '# cm^2 of N2');


end


% Define the concentration of NO2
% --------------------------------------------------------------------
if inputs.RT.modify_NO2==true

    % If true, modify the mixing ratio of molecular nitrogen
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s %s \n\n';
    fprintf(fileID, formatSpec,'mol_modify NO2 ', inputs.RT.NO2_mixing_ratio, 'CM_2', ' ', '# cm^2 of NO2');


end


% Define the concentration of molecular oxygen
% --------------------------------------------------------------------
if inputs.RT.modify_O2==true

    % If true, modify the mixing ratio of molecular oxygen
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'mixing_ratio O2 ', inputs.RT.O2_mixing_ratio, ' ', '# cm^2 of O2');


end


% Define the concentration of molecular oxygen
% --------------------------------------------------------------------
if inputs.RT.modify_O3==true

    % If true, modify the mixing ratio of ozone
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s %s \n\n';
    fprintf(fileID, formatSpec,'mol_modify O3 ', inputs.RT.O3_mixing_ratio, 'DU', ' ', '# Dobson Units of O3');


end



% Define the Aerosol Layer properties, if you want a cloud in your model
% --------------------------------------------------------------------
if inputs.RT.yesAerosols==true

    % Turn on default aersol layer, which occupies lower 2km of model
    % --------------------------------------------------------------
    formatSpec = '%s %5s %s \n';
    fprintf(fileID, formatSpec,'aerosol_default', ' ', '# turn on Shettle (1989) boundary layer aerosols');


    % Specify the Aerosl type
    % 1=rural aersols,  4=maritime aersols,  5=Urban aerosols,
    % 6=Tropospheric aerosols
    % ------------------------------------------------
    formatSpec = '%s %u %5s %s \n';
    fprintf(fileID, formatSpec,'aerosol_haze', inputs.RT.aerosol_type, ' ', '# Aerosol type');


    % Define aerosol layer optical depth
    % ----------------------------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'aerosol_modify tau set', inputs.RT.aerosol_opticalDepth, ' ',...
        '# Optical Depth of aerosol layer');

end



% Define the atm grid
% ------------------------------------------------
if isfield(inputs.RT, 'define_atm_grid') && inputs.RT.define_atm_grid==true
    % Define a custom atm grid
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'atm_z_grid', num2str(inputs.RT.atm_z_grid), ' ', '# Set vertical resolution of atm grid');

end



% Define the sensor altitude
% ------------------------------------------------
if ischar(inputs.RT.sensor_altitude)==true
    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'zout', inputs.RT.sensor_altitude, ' ', '# Sensor Altitude');

elseif isnumeric(inputs.RT.sensor_altitude)==true
    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'zout', num2str(inputs.RT.sensor_altitude), ' ', '# Sensor Altitude');

elseif iscell(inputs.RT.sensor_altitude)==true

    error([newline, 'Dont know how to deal with stings and numbers', newline])

    % find the string?
    idx_string = isstring(inputs.RT.sensor_altitude);
    idx_numeric = isnumeric(inputs.RT.sensor_altitude);
    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = ['%s %',num2str(length(idx_numeric)),'f %',num2str(length(idx_string)),'s %s \n'];
    fprintf(fileID, formatSpec,'zout', inputs.RT.sensor_altitude, ' ', '# Sensor Altitude');


end




% Define the solar zenith angle
% ------------------------------------------------
formatSpec = '%s %f %5s %s \n';
fprintf(fileID, formatSpec,'sza', inputs.RT.sza, ' ', '# Solar zenith angle');

% Define the solar azimuth angle
% -------------------------------------------------------
formatSpec = '%s %f %5s %s \n';
fprintf(fileID, formatSpec,'phi0', inputs.RT.phi0, ' ', '# Solar azimuth angle');

% Define the cosine of the zenith viewing angle
% ------------------------------------------------
formatSpec = '%s %f %5s %s \n';
fprintf(fileID, formatSpec,'umu', round(cosd(inputs.RT.vza),4), ' ', '# Cosine of the zenith viewing angle');

% Define the azimuth viewing angle
% ------------------------------------------------
formatSpec = '%s %f %5s %s \n\n';
fprintf(fileID, formatSpec,'phi', inputs.RT.vaz, ' ', '# Azimuthal viewing angle');



if inputs.RT.no_molecular_abs==true
    % turn off molecular absorption
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'no_absorption', 'mol', ' ', '# Turn off molecular absorption');
end


if inputs.RT.no_scattering_mol==true
    % Turn off molecular scattering
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'no_scattering', 'mol', ' ', '# Turn off molecular scattering');
end


if inputs.RT.no_scattering_aer==true
    % Turn off aerosol scattering
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'no_scattering', 'aer', ' ', '# Turn off aerosol scattering');
end


if inputs.RT.compute_reflectivity_uvSpec==true
    % Set the output quantity to be reflectivity
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'output_quantity', 'reflectivity', ' ', '# Output is reflectance');
end


%     % Set the outputs
%     % ------------------------------------------------
%     formatSpec = '%s %s %5s %s \n\n';
%     fprintf(fileID, formatSpec,'output_user', 'lambda edir edn eup uavgdir uavgdn uavgup uu', ' ', '# Output quantities');





% Set the error message to quiet of verbose
% ------------------------------------------------
formatSpec = '%s';
fprintf(fileID, formatSpec, inputs.RT.errMsg);


% Close the file!
fclose(fileID);
% ----------------------------------------------------
% ----------------------------------------------------






end