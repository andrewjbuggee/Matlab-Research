%% Define EMIT spectral response functions using the EMIT specified FWHM

% By Andrew John Buggee


function spec_response = create_EMIT_specResponse(emit, inputs)



spec_response.value = cell(length(emit.radiance.wavelength), 1);
spec_response.wavelength = cell(length(emit.radiance.wavelength), 1);

for ww = 1:length(emit.radiance.wavelength)

    % first we will create and store the spectral response function from
    % the full-wdith-half-max provided for each spectral channel
    % the spectral response function is a gaussian function
    % the emit wavelength vector is the center wavelength
    
    % set the center wavelength as the value defined by the EMIT wavelength
    % grid
    lambda_center = emit.radiance.wavelength(ww);      % nm

    % compute the standard deviation from the FWHM
    sigma = emit.radiance.fwhm(ww)/(2*sqrt(2*log(2)));      % std

    % create a wavelength vector
    if inputs.RT.source_file_resolution==0.1
        
        wl = round(lambda_center-(1.5*emit.radiance.fwhm(ww)), 1):...
            0.1:round(lambda_center+(1.5*emit.radiance.fwhm(ww)), 1);

    elseif inputs.RT.source_file_resolution==1
        
        wl = round(lambda_center-(1.5*emit.radiance.fwhm(ww))):...
            round(lambda_center+(1.5*emit.radiance.fwhm(ww)));

    end

    % compute the gaussian spectral response function
    spec_response.value{ww} = pdf('Normal', wl, lambda_center, sigma)';

    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    spec_response.wavelength{ww} = wl;

end





end