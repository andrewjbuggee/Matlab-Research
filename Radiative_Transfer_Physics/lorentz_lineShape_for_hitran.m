%% Compute the Pressure broadened line shape (Lorentz line shape)

% In the lower atmosphere, pressure broadening dominates


% ----------- INPUTS -----------------

% (1) hitarn_lines - the hitran output as define by the function
% importhitran()

% (2) wavelength_boundaries - (nm) the staring and ending wavelength as a 2
% element vector. 

% (3) T - (K) Temperature of the gas

% (4) P - (atm) is the total atmospheric pressure in atmospheres

% (5) P_self (atm) is the partial pressure in atmospheres corresponding to the molecule being calculated


% ----------- OUTPUTS -----------------

% (1) f_pressure - a structure containing the independent variable, which
% is the wavenumber encompassing the pressure broadened line, the
% dependent variable, the Lorentz line shape, and the
% half-width-at-half-max (HWHM)


% By Andrew John Buggee

%%

function f_pressure = lorentz_lineShape_for_hitran(hitran_lines, wavelength_boundaries, T, P, P_self)


% Specify the wavelength index using the range of interest
wl_index = hitran_lines.transitionWavenumber>=(10^4/(wavelength_boundaries(end)/1e3)) &...
    hitran_lines.transitionWavenumber<=(10^4/(wavelength_boundaries(1)/1e3));

% store variables you'll need in the for loop
transition_wavenumbers = hitran_lines.transitionWavenumber(wl_index);
pressure_shift = hitran_lines.pressureShift(wl_index);


% we have to compute a line shape for each spectral transition
lineShape_length = 100;

% define the number of HWHM for which we will compute the lineshape
num_hwhm = 11;


% define the hitran reference temperature
T_ref = 296;        % K

% compute the doppler profile half-width at half-max
% use SI units for everything, except leave the spectral dimension in
% wavenumbers with units of inverse centimeters
f_pressure.hwhm = (T_ref/T).^hitran_lines.temperatureDependence(wl_index) .* (hitran_lines.airBroadenedWidth(wl_index) .* (P - P_self) +...
    hitran_lines.selfBroadenedWidth(wl_index) .* P_self);


% compute the pressure broadened line shape for each line transition
f_pressure.shape = zeros(length(transition_wavenumbers), lineShape_length);
f_pressure.wavenum = zeros(length(transition_wavenumbers), lineShape_length);

for vv = 1:length(transition_wavenumbers)

    % Lines are narrow! Some less than a wavenumber
    f_pressure.wavenum(vv,:) = linspace(transition_wavenumbers(vv) - num_hwhm*f_pressure.hwhm(vv),...
        transition_wavenumbers(vv) + num_hwhm*f_pressure.hwhm(vv), lineShape_length);


    f_pressure.shape(vv,:) = f_pressure.hwhm(vv)./pi ./...
        (f_pressure.hwhm(vv).^2 + (f_pressure.wavenum(vv,:) - (transition_wavenumbers(vv) +...
        pressure_shift(vv) .* P)).^2 );

end



end