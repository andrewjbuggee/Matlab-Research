%% Define EMIT spectral response functions using the EMIT specified FWHM

% By Andrew John Buggee


function [inputs, spec_response] = create_EMIT_specResponse(emit, inputs)


% -------------------------------------------------------------------
% --- Compute using the EMIT center wavelengths and an arbitray width ---
% -------------------------------------------------------------------

% spec_response.value = cell(length(emit.radiance.wavelength), 1);
% spec_response.wavelength = cell(length(emit.radiance.wavelength), 1);
%
% for ww = 1:length(emit.radiance.wavelength)
%
%     % first we will create and store the spectral response function from
%     % the full-wdith-half-max provided for each spectral channel
%     % the spectral response function is a gaussian function
%     % the emit wavelength vector is the center wavelength
%
%     % set the center wavelength as the value defined by the EMIT wavelength
%     % grid
%     lambda_center = emit.radiance.wavelength(ww);      % nm
%
%     % compute the standard deviation from the FWHM
%     sigma = emit.radiance.fwhm(ww)/(2*sqrt(2*log(2)));      % std
%
%     % create a wavelength vector using the source file resolution
%     source_file_resolution = inputs.RT.source_file_resolution;      % nm
%
%     % The wavelength vector for libRadTran is simply the lower and upper
%     % bounds
%     spec_response.wavelength{ww}  = round(lambda_center-(1.5*emit.radiance.fwhm(ww)), 1):...
%          source_file_resolution:round(lambda_center+(1.5*emit.radiance.fwhm(ww)), 1);
%
%     % compute the gaussian spectral response function
%     spec_response.value{ww} = pdf('Normal', spec_response.wavelength{ww}, lambda_center, sigma)';
%
%
% end



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
[~, wavelength_source] = read_solar_flux_file([300, 2600], inputs.RT.source_file);

spec_response.value = zeros(length(emit.radiance.wavelength), length_spec_response_func+1);
spec_response.wavelength = zeros(length(emit.radiance.wavelength), length_spec_response_func+1);

for ww = 1:length(emit.radiance.wavelength)

    % In libRadtran, everything is computed on the source wavelength grid.
    % So if the source is defined at integer values of wavelength (100,
    % 101, ...), then the computations must be done at those wavelengths

    % Find solar flux wavelength closest to each center wavelength of the
    % EMIT spectrometer
    [~, center_wavelength_idx] = min(abs(wavelength_source - emit.radiance.wavelength(ww)));

    % set the center wavelength as the value defined by the wavelength
    % closest to the EMIT center wavelength
    lambda_center = wavelength_source(center_wavelength_idx);      % nm

    % compute the standard deviation from the FWHM
    sigma = emit.radiance.fwhm(ww)/(2*sqrt(2*log(2)));      % std

    % define the wavelength grid for each spectral channel
    spec_response.wavelength(ww, :) = wavelength_source(center_wavelength_idx - length_spec_response_func/2:...
        center_wavelength_idx + length_spec_response_func/2);

    % compute the gaussian spectral response function
    spec_response.value(ww, :) = pdf('Normal', spec_response.wavelength(ww, :), lambda_center, sigma);




end


% Now, define the wavelength range for each channel being modeled for EMIT
% now define the wavelength range of each spectral channel
inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);

for ww = 1:numel(inputs.bands2run)

    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    inputs.RT.wavelengths2run(ww,:) = [spec_response.wavelength(inputs.bands2run(ww), 1),...
        spec_response.wavelength(inputs.bands2run(ww), end)];

end




end