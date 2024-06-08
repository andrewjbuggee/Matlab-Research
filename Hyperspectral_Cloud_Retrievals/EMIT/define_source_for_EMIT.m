%% Compute the at-sensor TOA solar flux as measured by EMIT

% We compute a vector that is the same length as the EMIT measurements in spectral space.  
% This is accomplished by integrating the TOA solar flux over each EMIT spectral channel
% using the spectral response function

% ----- OUTPUTS -----

% (1) source.flux - W/nm/m^2
% (2) source.wavelength - nm

% By Andrew John Buggee



%%

function inputs = define_source_for_EMIT(inputs, emit)


% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------
% resolution should match the value listed in the file name
source_file_resolution = inputs.RT.source_file_resolution;
source_file = inputs.RT.source_file;

% 
% if source_file_resolution==0.1
%     
%     inputs.RT.source_file = 'kurudz_0.1nm.dat';
% 
% 
% elseif source_file_resolution==1
% 
%     inputs.RT.source_file = 'kurudz_1.0nm.dat';
% 
% 
% elseif source_file_resolution==0.025
% 
%     inputs.RT.source_file = 'hybrid_reference_spectrum_p025nm_resolution_c2022-11-30_with_unc.nc';
%     
% 
% elseif source_file_resolution==0.005
% 
%     inputs.RT.source_file = 'hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc';
% 
% else
% 
%     error([newline, 'What wavelength resolution do you wish to use?', newline])
% 
% end



% read the source file
% [source.flux, source.wavelength] = read_solar_flux_file([emit.spec_response.wavelength{1}(1), ...
%     emit.spec_response.wavelength{end}(end)], inputs.RT.source_file);        % (W/nm/m^2) 

[source.flux, source.wavelength] = read_solar_flux_file([300, 2600], source_file);        % (W/nm/m^2) 

% Trim the solar flux data so it uses the same wavelength gird as EMIT
% Do this by integrating the spectral response function of each
% spectral channel with the solar flux file


% 
% ------------------------------------------------------------------------
% Integrate the spectral response functions over the EMIT center wavelength
% ------------------------------------------------------------------------
% 
% inputs.source.flux = zeros(length(emit.radiance.wavelength), 1);
% inputs.source.wavelength = zeros(length(emit.radiance.wavelength), 1);
% 
% for ww = 1:length(emit.radiance.wavelength)
% 
%     %disp(num2str(ww))
%     % Find solar flux wavelength  closest to each center wavelength of the
%     % EMIT spectrometer
%     [~, center_wavelength_idx] = min(abs(source.wavelength - emit.radiance.wavelength(ww)));
% 
%     % set the center wavelength as the value defined by the wavelength
%     % closest to the EMIT center wavelength
%     lambda_center = source.wavelength(center_wavelength_idx);      % nm
% 
%     % compute the standard deviation from the FWHM
%     sigma = emit.radiance.fwhm(ww)/(2*sqrt(2*log(2)));      % std
% 
%     % define the wavelength grid for each spectral channel
%     wl = source.wavelength(center_wavelength_idx - 600 : center_wavelength_idx + 600);
% 
%     % compute the gaussian spectral response function
%     spec_response = pdf('Normal', wl, lambda_center, sigma);
% 
%     % define an index in order to use the appropriate soalr flux values
%     wavelength_idx = source.wavelength >= wl(1) & source.wavelength <= wl(end);
% 
%     % integrate the source flux with the spectral response function
%     inputs.source.flux(ww) = trapz(wl, spec_response .* source.flux(wavelength_idx));           % W/nm/m^2
% 
%     inputs.source.wavelength(ww) = lambda_center;    % nm
% 
% 
% %     wavelength_idx = source.wavelength >= emit.spec_response.wavelength{ww}(1) &...
% %         source.wavelength <= emit.spec_response.wavelength{ww}(end);
% % 
% %     inputs.source.flux(ww) = trapz(source.wavelength(wavelength_idx),...
% %         emit.spec_response.value{ww} .* source.flux(wavelength_idx));           % W/nm/m^2
% % 
% %     inputs.source.wavelength(ww) = source.wavelength(center_wavelength_idx);    % nm
% 
% end


% ------------------------------------------------------------------------
% Convolve the spectral response functions over the solar spectral flux
% ------------------------------------------------------------------------


% pick a FWHM from the middle of the spectrum
% compute the standard deviation from the FWHM
sigma = emit.radiance.fwhm(136)/(2*sqrt(2*log(2)));      % std

% compute the gaussian spectral response function
spec_response = pdf('Normal', source.wavelength, median(source.wavelength), sigma);

% Convolve each point of the solar flux grid with the spectral response
% function

convolved_solar_flux = conv(source.flux, spec_response, 'same') .* source_file_resolution;


% interpolate the convolved solar flux to get the values on the EMIT
% wavelength grid
inputs.source.flux = interp1(source.wavelength, convolved_solar_flux, emit.radiance.wavelength);    % W/m^2/nm
inputs.source.wavelength = emit.radiance.wavelength;        %nm







% % define the length of the solar flux wavelength grid
% num_wl = length(source.wavelength);
% 
% % define the boundaries of the solar_flux wavelength grid
% wl_min = source.wavelength(1);      % nm
% wl_max = source.wavelength(end);    % nm
%
% % define the new spectrum, the solar flux convolved with the spectral
% % response
% convolved_solar_flux = zeros(length(source.wavelength), 1);
% 
% % step through each wavelength point on the solar flux grid
% for ww = 1:length(source.wavelength)
%     
%     % define the spectral response function so that it is twice the length of
%     % the solar flux wavelength grid
%     wl_spec_response_grid = ((source.wavelength(ww) - num_wl*source_file_resolution):...
%         source_file_resolution:(source.wavelength(ww) + num_wl*source_file_resolution))';   % nm
% 
%     % append zeros on either side of the solar flux so it is the same length as
%     % the wl_spec_response_grid
%     % first find the wavelength values outside the solar flux wavelength
%     % grid
% 
%     index_outOfBounds = wl_spec_response_grid>=wl_min & wl_spec_response_grid<=wl_max;
%     solar_flux = zeros(length(index_outOfBounds), 1);
%     solar_flux(index_outOfBounds) = source.flux;
% 
%     % create the spectral response function
%     % set the center wavelength as the value defined by the EMIT wavelength
%     % grid
%     lambda_center = source.wavelength(ww);
% 
%     % compute the gaussian spectral response function
%     spec_response = pdf('Normal', source.wavelength, median(source.wavelength), sigma);
% 
%     % Convolve each point of the solar flux grid with the spectral response
%     % function
%     
%     convolved_solar_flux = conv(source.flux, spec_response, 'same') .* source_file_resolution;
% 
%     inputs.source.flux(ww) = trapz(source.wavelength(wavelength_idx),...
%         emit.spec_response.value{ww} .* source.flux(wavelength_idx));           % W/nm/m^2
% 
%     inputs.source.wavelength(ww) = source.wavelength(center_wavelength_idx);    % nm
% 
% end



end
