% ----- Planks Function -----

% this function takes in a vector of wavelength values, and a vector of
% temperature values and computes the plank function, which defines the
% amount of radiance 

% ----- inputs -----

% (1) - independentVariable - a row vector defining the independent
% variable used for the planck function. The choice is arbitrary. One can
% define the planck function so that it varies with wavelength, frequency,
% or wavenumber. Though I don't know why you'd want to use wavenumbers. 

% (2) - temperature - a row vector. For each temperature, we create a
% row vector of radiance across the independent variable defined by the user.
% units are in kelvin. The temperature of the photosphere of the sun is
% 5778 K

% (3) - indVar_units - a string defining the units to use for both the
% input and the output. If you are using wavelength as the independent
% variable, you can use either:
%       (a) 'microns'
%       (b) 'nanometers'

% If you're using frequency as the independent variable, you have only 
% one option for units 
%       (a) 'Hz'

% If you're using wavenumber as the independent variable, you have only 
% one option for units 
%       (a) 'cm^-1'


% ----- outputs -----

% The units of the output depends on the defined independent variable units

% (1) radiance - units are in watts/m^2/nm/sr
% or             units are in watts/m^2/micron/sr
% or             units are in watts/m^2/Hz/sr
% or             units are in watts/m^2/cm^{-1}/sr
                    


% By Andrew J. Buggee
%%

function [radiance] = plancks_function(independent_variables,temperatures,indVar_units)

% define constants

h = 6.626e-34; % Joules * sec - planks constant
c = 299792458; % meters/s - speed of light
k = 1.3806e-23; % J/Kelvin - Boltzmann constant


% output a matrix with number of rows equal to the number of temperatures,
% and the number of columns equal to the number of spectral bins
radiance = zeros(length(temperatures),length(independent_variables));

if strcmp(indVar_units,'microns')
    
    % convert wavelenghts to meters
    independent_variables = independent_variables .* 1e-6; % meters
    
    for ii = 1:length(temperatures)
        
        radiance(ii,:) = 2*h*c^2./(independent_variables.^5 .*(exp(h*c./(independent_variables*k*temperatures(ii))) - 1));
        
    end
    
    % radiance now has units of W/m^3/sr. We need to convert this to be
    % W/m^2/micron/sr
    
    radiance = radiance.* 1/1e6; % - W/m^2/micron/sr
    
    
elseif strcmp(indVar_units,'nanometers')
    
    % convert wavelenghts to meters
    independent_variables = independent_variables .* 1/1e9; % meters
    
    for ii = 1:length(temperatures)
        
        radiance(ii,:) = 2*h*c^2./(independent_variables.^5 .*(exp(h*c./(independent_variables*k*temperatures(ii))) - 1));
        
    end
    
    % radiance now has units of W/m^3/sr. We need to convert this to be
    % W/m^2/nm/sr
    
    radiance = radiance.* 1/1e9; % - W/m^2/nm/sr


elseif strcmp(indVar_units,'Hz')
    
    % the independent variable is defined in units of Hertz, or inverse
    % seconds

    % This requires a different definition of the planck function, where
    % one must transform the wavelength definition using the jacobian,
    % dL/dv - the change in wavelength versus the change in frequency
    
    for ii = 1:length(temperatures)
        
        % radiance has units of W/m^2/sr/Hz. 
        radiance(ii,:) = 2*h*independent_variables.^3./(c^2 .*(exp(h*independent_variables./(k*temperatures(ii))) - 1));  % - W/m^2/Hz/sr
        
    end
    


elseif strcmp(indVar_units,'cm^-1')
    
    % the independent variable is defined in units of inverse cm, also
    % known as wavenumbers
    
    % convert inverse cm to inverse m
    independent_variables = independent_variables * 100;            % m^{-1}

    % This requires a different definition of the planck function, where
    % one must transform the wavelength definition using the jacobian,
    % dL/dv' - the change in wavelength versus the change in wavenumber
    
    for ii = 1:length(temperatures)
        
        radiance(ii,:) = 2*h*c^2*independent_variables.^3./(exp(h*c*independent_variables./(k*temperatures(ii))) - 1);
        
    end
    
    % radiance now has units of W/m^2/m^{-1}/sr. We need to convert this to be
    % W/m^2/cm^{-1}/sr
    
    radiance = radiance * 100; % - W/m^2/cm^{-1}/sr

    
else
    
    error([newline,'I dont recognize your input for units. Acceptable inputs are: "microns"',...
        ' or "nanometers".'])
    
end




end

