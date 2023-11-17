%% This function computes brightness temperturature

% Assuming an emissivity of 1, this functions computes the temperature a
% black body would need to have in order to produce some spectral output,
% L. The derivative of the Planck function with respect to temperature is
% monotonic with semi-infinite bounds. Thus there is always a unique
% brightness temperature that satisfies our spectral measurements.


% ----- inputs -----

% (1) radiance_measurement - This is the radiance used to compute
% brightness temperature. The input can be a value at a monochromatic
% wavelength, or a radiance measurement over some spectral region. The
% latter is how real instruments before brightness temperature
% calculations, since all real instruments have to integrate over some
% finite spectral window. Only a single value is required.

% (2) - independentVariable - Define the independent variable used to make 
% the spectral measurement. If using a monochromatic measurement, only a
% single value is required. If using a real measurement over some spectral
% window, from the spectral bounds of the measurement. The choice is 
% arbitrary. One can define the planck function so that it varies with 
% wavelength, frequency, or wavenumber. Though I don't know why you'd
% want to use wavenumbers. 


% (3) - indVar_units - a string defining the units to use for both the
% input and the output. If you are using wavelength as the independent
% variable, you can use either:
%       (a) 'microns'
%       (b) 'nanometers'

    % If you're using frequency as the independent variable, you have 
    % only one option for units 
    %       (a) 'Hz'

    % If you're using wavenumber as the independent variable, you have 
    % only one option for units 
    %       (a) 'cm^-1'


% ----- outputs -----


% (1) Brightness Temperature (Tb) - units in Kelvin (K)


% By Andrew John Buggee
%%

function Tb = brightness_temperature(radiance_measurement, independentVariable, indVar_units)

con = physical_constants();

if strcmp(indVar_units, 'microns')==true


    % Check to see if the measurement is monochromatic or integrated over
    % some spectral region

    if length(independentVariable)==1

        spec_var = independentVariable./1e6;                      % meters - converted from microns
        L = radiance_measurement.*1e6;                            % W/m^2/sr/meter - converted to a spectral dependence using meters

    
        Tb = con.h * con.c./(spec_var * con.k_B) * 1/(log(2*con.h * con.c^2/(spec_var.^5 .*L) +1));

    
    elseif length(independentVariable)==2

        % Then our input simulates a real measurement, and we have to
        % integrate over the spectral window provided.
        % We cannot directly solve for a temperature if the measurement is
        % integrated over some spectral window. Instead we have to find the
        % temperature that provides a measurement close to the input

        % The resolution will be a 1 degree K
        temperature = 1:10000;          % K - temperatures to test

        % create a smooth independent variable to integrate over
        spectral_variable = linspace(independentVariable(1), independentVariable(2), 250);
        
        planck_measurements = trapz(spectral_variable, plancks_function(spectral_variable,temperature,indVar_units),2);

        % find the minimum absolute distance
        [~, min_dist] = min(abs(planck_measurements - radiance_measurement));

        % find the corresponding temperature
        Tb = temperature(min_dist);             % K
        

    else

        error([newline,'The spectral variable input can only be of length 1 or 2.', newline])

    end




elseif strcmp(indVar_units, 'nanometers')


     % Check to see if the measurement is monochromatic or integrated over
    % some spectral region

    if length(independentVariable)==1

        spec_var = independentVariable./1e9;                      % meters - converted from microns
        L = radiance_measurement.*1e9;                            % W/m^2/sr/meter - converted to a spectral dependence using meters

    
        Tb = con.h * con.c./(spec_var * con.k_B) * 1/(log(2*con.h * con.c^2/(spec_var.^5 .*L) +1));     % K

    
    elseif length(independentVariable)==2

        % Then our input simulates a real measurement, and we have to
        % integrate over the spectral window provided.
        % We cannot directly solve for a temperature if the measurement is
        % integrated over some spectral window. Instead we have to find the
        % temperature that provides a measurement close to the input

        % The resolution will be a 1 degree K
        temperature = 1:10000;          % K - temperatures to test

        % create a smooth independent variable to integrate over
        spectral_variable = linspace(independentVariable(1), independentVariable(2), 250);
        
        planck_measurements = trapz(spectral_variable, plancks_function(spectral_variable,temperature,indVar_units),2);

        % find the minimum absolute distance
        [~, min_dist] = min(abs(planck_measurements - radiance_measurement));

        % find the corresponding temperature
        Tb = temperature(min_dist);             % K
        

    else

        error([newline,'The spectral variable input can only be of length 1 or 2.', newline])

    end




else

    error([newline,'I dont recognize the spectral units of your input.',newline])

end



end