%% Write the INP files used by libRadtran to solve the radiative transfer equation

% By Andrew John Buggee
%%

function [] = write_INP_file(INP_folderpath, libRadtran_data_path, inputFileName, inputs, wavelengths, wc_filename)



% ----------------- ******************** -------------------
% ------------------ Write the INP File --------------------
% ----------------- ******************** -------------------

% Open the old file for writing
fileID = fopen([INP_folderpath, inputFileName], 'w');

% Define which RTE solver to use
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n';
fprintf(fileID, formatSpec,'rte_solver', inputs.RT.rte_solver,' ', '# Radiative transfer equation solver');


% Define the number of streams to keep track of when solving the equation
% of radiative transfer
% ------------------------------------------------
formatSpec = '%s %u %5s %s \n\n';
fprintf(fileID, formatSpec,'number_of_streams', inputs.RT.num_streams,' ', '# Number of streams');


% Use phase function correction?
% ------------------------------------------------
if inputs.RT.use_nakajima_phaseCorrection==true
    % define the pahse correction to be true
    % ------------------------------------------------
    formatSpec = '%s %5s %s \n\n';
    fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
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



% Define the surface albedo
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n\n';
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



% Define the column water vapor amount
% --------------------------------------------------------------------
if inputs.RT.modify_waterVapor==true

    % If true, modify the amount of column water vapor
    % --------------------------------------------------------------
    formatSpec = '%s %f %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'mol_modify H2O ', inputs.RT.waterVapor_column, ' MM', ' ', '# Column water vapor amount');


end


% Define the concentration of carbon dioxide
% --------------------------------------------------------------------
if inputs.RT.modify_CO2==true

    % If true, modify the mixing ratio of carbon dioxide
    % --------------------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'mixing_ratio CO2 ', inputs.RT.CO2_mixing_ratio, ' ', '# ppm of CO2');


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




% Define the sensor altitude
% ------------------------------------------------
formatSpec = '%s %s %5s %s \n';
fprintf(fileID, formatSpec,'zout', inputs.RT.sensor_altitude, ' ', '# Sensor Altitude');

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