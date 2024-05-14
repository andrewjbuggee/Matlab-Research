%% Compute the Pressure broadened line shape (Lorentz line shape)

% In the lower atmosphere, pressure broadening dominates


% ----------- INPUTS -----------------

% (1) hitarn_lines - the hitran output as define by the function
% importhitran()

% (2) T - (K) Temperature of the gas

% (3) P - (atm) is the total atmospheric pressure in atmospheres

% (4) P_self (atm) is the partial pressure in atmospheres corresponding to the molecule being calculated


% ----------- OUTPUTS -----------------

% (1) f_pressure - a structure containing the independent variable, which
% is the wavenumber encompassing the pressure broadened line, and the
% dependent variable, the Lorentz line shape.


% By Andrew John Buggee

%%

function f_pressure = lorentz_lineShape_for_hitran(hitran_lines, T, P, P_self)

% we have to compute a line shape for each spectral transition
lineShape_length = 100;

% define the number of HWHM for which we will compute the lineshape
num_hwhm = 11;


% define the hitran reference temperature
T_ref = 296;        % K

% compute the doppler profile half-width at half-max
% use SI units for everything, except leave the spectral dimension in
% wavenumbers with units of inverse centimeters
pressure_hwhm = (T_ref/T).^hitran_lines.temperatureDependence .* hitran_lines.airBroadenedWidth .* (P - P_self) +...
    hitran_lines.selfBroadenedWidth .* P_self;


% compute the pressure broadened line shape for each line transition
f_pressure.shape = zeros(length(hitran_lines.transitionWavenumber), lineShape_length);
f_pressure.wavenum = zeros(length(hitran_lines.transitionWavenumber), lineShape_length);

for vv = 1:length(hitran_lines.transitionWavenumber)

    % Lines are narrow! Some less than a wavenumber
    f_pressure.wavenum(vv,:) = linspace(hitran_lines.transitionWavenumber(vv) - num_hwhm*pressure_hwhm(vv),...
        hitran_lines.transitionWavenumber(vv) + num_hwhm*pressure_hwhm(vv), lineShape_length);


    f_pressure.shape(vv,:) = pressure_hwhm(vv)./pi ./...
        (pressure_hwhm(vv).^2 + (f_pressure.wavenum(vv,:) - (hitran_lines.transitionWavenumber(vv) +...
        hitran_lines.pressureShift(vv) .* P)).^2 );

end



end