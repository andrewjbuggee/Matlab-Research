%% Compute the at-sensor TOA solar flux as measured by EMIT

% We compute a vector that is the same length as the EMIT measurements in spectral space.  
% This is accomplished by integrating the TOA solar flux over each EMIT spectral channel
% using the spectral response function

% By Andrew John Buggee



%%

function inputs = define_source_for_EMIT(inputs, emit)


% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------
% resolution should match the value listed in the file name

if inputs.RT.source_file_resolution==0.1
    
    inputs.RT.source_file = '../data/solar_flux/kurudz_0.1nm.dat';

elseif source_file_resolution==1

    inputs.RT.source_file = '../data/solar_flux/kurudz_1.0nm.dat';

end



% read the source file
[source.flux, source.wavelength] = read_solar_flux_file([emit.spec_response.wavelength{1}(1), ...
    emit.spec_response.wavelength{end}(end)], inputs.RT.source_file(20:end));        % (W/nm/m^2) 

% Trim the solar flux data so it is the same length as the EMIT wavelength
% vector. Do this by integrating the spectral response function of each
% spectral channel with the solar flux file

inputs.source.flux = zeros(length(emit.radiance.wavelength), 1);
inputs.source.wavelength = zeros(length(emit.radiance.wavelength), 1);


for ww = 1:length(emit.radiance.wavelength)
    
    % Find solar flux wavelength closest to each center wavelength of the
    % EMIT spectrometer
    [~, center_wavelength_idx] = min(abs(source.wavelength - emit.radiance.wavelength(ww)));

    wavelength_idx = source.wavelength >= emit.spec_response.wavelength{ww}(1) &...
        source.wavelength <= emit.spec_response.wavelength{ww}(end);

    inputs.source.flux(ww) = trapz(source.wavelength(wavelength_idx),...
        emit.spec_response.value{ww} .* source.flux(wavelength_idx));           % W/nm/sr

    inputs.source.wavelength(ww) = source.wavelength(center_wavelength_idx);    % nm

end





end
