%% Define arbitrary Gaussian spectral response functions 

% By Andrew John Buggee


function spec_response = create_gaussian_specResponse(center_wavelength, fwhm, inputs)




% -------------------------------------------------------------------
% --- Another way of solving this is keeping everything on the source
% function wavelength grid ---
% -------------------------------------------------------------------

% define the length of the wavelength grid
if inputs.RT.source_file_resolution==1
    length_spec_response_func = 100;

elseif  inputs.RT.source_file_resolution==0.1
    length_spec_response_func = 300;

elseif inputs.RT.source_file_resolution==0.025
    length_spec_response_func = 600;

elseif inputs.RT.source_file_resolution==0.005
    length_spec_response_func = 800;

end

% read in the desired source file wavelength grid
[~, wavelength_source] = read_solar_flux_file([300, 2600], inputs.RT.source_file);          % nm - 

spec_response.value = zeros(length(center_wavelength), length_spec_response_func+1);
spec_response.wavelength = zeros(length(center_wavelength), length_spec_response_func+1);

for ww = 1:length(center_wavelength)

    % Find solar flux wavelength closest to each center wavelength
    [~, center_wavelength_idx] = min(abs(wavelength_source - center_wavelength(ww)));

    % set the center wavelength as the value defined by the wavelength
    % closest to the center wavelength
    lambda_center = wavelength_source(center_wavelength_idx);      % nm

    % compute the standard deviation from the FWHM
    sigma = fwhm(ww)/(2*sqrt(2*log(2)));      % std

    % define the wavelength grid for each spectral channel
    spec_response.wavelength(ww, :) = wavelength_source(center_wavelength_idx - length_spec_response_func/2:...
        center_wavelength_idx + length_spec_response_func/2);

    % compute the gaussian spectral response function
    spec_response.value(ww, :) = pdf('Normal', spec_response.wavelength(ww, :), lambda_center, sigma);


end


end