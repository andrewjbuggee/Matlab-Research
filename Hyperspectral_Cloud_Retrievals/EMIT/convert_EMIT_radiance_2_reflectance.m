%% Convert Measured at-sensor EMIT radiance to TOA reflectance


% *** Spectral response function not needed because it is already
% incorporated with the radiance measurements and the solar flux ***


% By Andrew John Buggee
%%

function emit = convert_EMIT_radiance_2_reflectance(emit, inputs)


% read in the source file
% The solar flux, as defined by define_source_for_emit.m, is the same
% length as the spectral dimension of the EMIT data cube.

% *** IMPORTANT *** the solar flux has been integrated over each EMIT
% spectral channel with the spectral response function
solar_flux = inputs.source.flux;     % W/m^2/nm

% Convert this into the same units used by EMIT (microW/cm^2/nm)
solar_flux = solar_flux .* 10^2;        % microW/cm^2/nm


% Compute the reflectance
% if there is only 1 pixel, we don't need a for loop

if length(emit.obs.solar.zenith)==1

    % Compute EMIT measurements to reflectance according to eq. 1 from
    % Thompson et al. 2016
    emit.reflectance = pi*emit.radiance.measurements ./...
        (cosd(emit.obs.solar.zenith) * solar_flux);             % 1/sr


else

    emit.reflectance = zeros(length(emit.radiance.wavelength), length(emit.obs.solar.zenith));

    % Step through each pixel and compute the reflectance
    for nn = 1:length(emit.obs.solar.zenith)

        % Compute EMIT measurements to reflectance according to eq. 1 from
        % Thompson et al. 2016
        emit.reflectance(:, nn) = pi*emit.radiance.measurements(:,nn) ./...
            (cosd(emit.obs.solar.zenith(nn)) * solar_flux);             % 1/sr

    end



end







end