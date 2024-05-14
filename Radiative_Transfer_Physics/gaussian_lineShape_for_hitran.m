%% Compute the Doppler broadened line shape (gaussian line shape)

% Doppler broadening tends to dominate in the upper atmosphere where
% pressure is low

% ----------- INPUTS -----------------

% (1) hitarn_lines - the hitran output as define by the function
% importhitran()

% (2) T - (K) Temperature of the gas


% ----------- OUTPUTS -----------------

% (1) f_doppler - a structure containing the independent variable, which
% is the wavenumber encompassing the doppler broadened line, and the
% dependent variable, the Gaussian line shape.


% By Andrew John Buggee

%%

function f_doppler = gaussian_lineShape_for_hitran(hitran_lines, T)

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
doppler_hwhm = hitran_lines.transitionWavenumber./con.c .* sqrt((2*con.N_A * con.k_B * T * log(2))/molar_mass);


% compute the Doppler line shape for each line transition
f_doppler.shape = zeros(length(hitran_lines.transitionWavenumber), lineShape_length);
f_doppler.wavenum = zeros(length(hitran_lines.transitionWavenumber), lineShape_length);

for vv = 1:length(hitran_lines.transitionWavenumber)

    % Lines are narrow! Some less than a wavenumber
    f_doppler.wavenum(vv,:) = linspace(hitran_lines.transitionWavenumber(vv) - num_hwhm*doppler_hwhm(vv),...
        hitran_lines.transitionWavenumber(vv) + num_hwhm*doppler_hwhm(vv), lineShape_length);

    f_doppler.shape(vv,:) = sqrt(log(2)./(pi*doppler_hwhm(vv).^2)) .* ...
        exp(- ((f_doppler.wavenum(vv,:) - hitran_lines.transitionWavenumber(vv)).^2 * log(2))./doppler_hwhm(vv)^2);

end



end