%% Compute the Doppler broadened line shape (gaussian line shape)

% Doppler broadening tends to dominate in the upper atmosphere where
% pressure is low

% ----------- INPUTS -----------------

% (1) hitarn_lines - the hitran output as define by the function
% importhitran()

% (2) wavelength_boundaries - (nm) the staring and ending wavelength as a 2
% element vector. 

% (3) T - (K) Temperature of the gas


% ----------- OUTPUTS -----------------

% (1) f_doppler - a structure containing the independent variable, which
% is the wavenumber encompassing the doppler broadened line, the
% dependent variable, the Gaussian line shape, and the
% hwafl-width-at-half-max (HWHM)


% By Andrew John Buggee

%%

function f_doppler = gaussian_lineShape_for_hitran(hitran_lines, wavelength_boundaries, T)


% Specify the wavelength index using the range of interest
wl_index = hitran_lines.transitionWavenumber>=(10^4/(wavelength_boundaries(end)/1e3)) &...
    hitran_lines.transitionWavenumber<=(10^4/(wavelength_boundaries(1)/1e3));

transition_wavenumbers = hitran_lines.transitionWavenumber(wl_index);



% we have to compute a line shape for each spectral transition
lineShape_length = 100;

% define the number of HWHM for which we will compute the lineshape
num_hwhm = 3;

% We need the molar mass of the isotopologue

molar_mass = read_isotopologue_molar_mass_hitran(hitran_lines);     % g/mol
% convert molar_mass to kg/mol
molar_mass = molar_mass/1000;                           % kg/mol

% load physical constants
con = physical_constants();

% compute the doppler profile half-width at half-max
% use SI units for everything, except leave the spectral dimension in
% wavenumbers with units of inverse centimeters
f_doppler.hwhm = transition_wavenumbers./con.c .* sqrt((2*con.N_A * con.k_B * T * log(2))/molar_mass);


% compute the Doppler line shape for each line transition



f_doppler.shape = zeros(length(transition_wavenumbers), lineShape_length);
f_doppler.wavenum = zeros(length(transition_wavenumbers), lineShape_length);

for vv = 1:length(transition_wavenumbers)

    % Lines are narrow! Some less than a wavenumber
    f_doppler.wavenum(vv,:) = linspace(transition_wavenumbers(vv) - num_hwhm*f_doppler.hwhm(vv),...
        transition_wavenumbers(vv) + num_hwhm*f_doppler.hwhm(vv), lineShape_length);

    f_doppler.shape(vv,:) = sqrt(log(2)./(pi*f_doppler.hwhm(vv).^2)) .* ...
        exp(- ((f_doppler.wavenum(vv,:) - transition_wavenumbers(vv)).^2 * log(2))./f_doppler.hwhm(vv)^2);

end



end