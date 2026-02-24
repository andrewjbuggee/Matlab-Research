%% ----- CALCULATE THE REFLECTANCE FUNCTION -----

% INPUTS:
%   (1) inputSettings

%   (2) dataStructure  - output results from uvspec

%   (3) spec_response - spectral response function - we need only the
%   weightings, not the wavelength data


% OUTPUTS:
%   (1) R - This is the reflectance over a finite bandwidth. R integrates
%   the monochromatic reflectance, R_lambda, over the wavelength range

%   (2) R_lambda - This is the reflectance for a monochromatic calculation.
%   If the equation of radiative transfer is solved for a single
%   wavelength, then we don't need to integrate. 


% By Andrew J. Buggee
%% ------ Read input settings and output data from uv_spec -----

function [R,R_lambda] = reflectanceFunction_ver3(ds, source, spec_response, sza, vza, vaz)


% Geometry values from input Settings -
mu = cosd(vza); % cosine of the viewing zenith angle





% radiative transfer solutions
wavelength = ds.wavelength;     % nm

% convert source to mW/(m^2 nm) -  source irradiance
source = source .* 1e3;            % - mW/(m^2 nm) 

% check to see if the source is a different length than the radiance vector
if length(source)~=length(ds.radiance.value)

    error([newline, 'The source file doesnt match the length of the computed radiance!', newline]);

end

% a few other constants calculated from the inputs
mu0 = cosd(sza); % cosine of the solar zenith angle
geomSets = length(mu)*length(vaz);


% make sure the spectral response is a column vector
spec_response = reshape(spec_response, [], 1);



%% ----- CALCULATE REFLECTANCE FUNCTION -----

% reflectance varies with lambda, tau_cloud, r_e, viewing angle, solar
% zenith angle, and the viewing azimuth angle. uvspec can only vary the
% viewing azimuth and viewing zenith angle for a given file. So we will
% create a 3D array that varies these three parameters per file

% ***** Radiance has units of mW/m^2/nm/sr *****

R_lambda = zeros(length(wavelength),geomSets); % phi (azimuth) changes first, then mu (cos(sza))
R = zeros(length(mu),length(vaz));
% if there is a single monochromatic wavelength then we don't need to
% integrate. We simply divide
if isscalar(wavelength)
    
    for ii = 1:geomSets
        
        R_lambda(ii) = pi*ds.radiance(ii).value./(mu0*source); % - 1/sr/nm - reflectance function for monochromatic calculation
        
    end
    
elseif length(wavelength)>1
    
    for ii = 1:geomSets
        
        % First calculate the reflectance function at each discrete
        % wavelength within the wavelength band
        R_lambda(:,ii) = pi*ds.radiance(ii).value./(mu0*source); % - 1/sr/nm - reflectance function for monochromatic calculation
        R(ii) = trapz(wavelength, R_lambda(:,ii).*spec_response.*source)./trapz(wavelength, spec_response.*source); % - 1/sr - reflectance function over a finite bandwidth
        
    end
    
else
    
    error('Something is wrong with the wavelength vector');
    
end









end