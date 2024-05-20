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

% (1) f_pressure - [(cm^(-1)]^-1 - a structure containing the independent variable, which
% is the wavenumber encompassing the pressure broadened line, the
% dependent variable, the Lorentz line shape, and the
% half-width-at-half-max (HWHM)


% By Andrew John Buggee

%%

function f_pressure = lorentz_lineShape_for_hitran(hitran_lines, wavelength_grid, T, P, P_self)


%% Convert the wavelength grid to a wavenumber grid


% ------------------------------------------------------------------------
% ************ DOESN'T PLACE LINE SHAPES ON THE SAME GRID ****************
% ------------------------------------------------------------------------


% Convert this linearly spaced wavelength grid into a wavenumber grid
% Make sure the wavelength vector is in microns
wavenumber_master_grid = 10^4 ./ (wavelength_grid./1e3);        % cm^(-1)

% Make sure wavenumber_master_grid is a row vector
if size(wavenumber_master_grid,1)==1 && size(wavenumber_master_grid,2)>1

    wavenumber_master_grid = wavenumber_master_grid';
end

% Find the line centers closest to each value of the wavenumber grid
[~, w_index] = min(abs(wavenumber_master_grid - hitran_lines.transitionWavenumber'), [], 2);

% Keep only the unique values
w_index = unique(w_index, 'stable');

% grab the the line strength centers
line_center = hitran_lines.transitionWavenumber(w_index);

%%
pressure_shift = hitran_lines.pressureShift(w_index);


% we have to compute a line shape for each spectral transition
lineShape_length = 100;

% define the number of HWHM for which we will compute the lineshape
num_hwhm = 11;


% define the hitran reference temperature
T_ref = 296;        % K

% compute the pressure broadened half-width at half-max
f_pressure.hwhm = (T_ref/T).^hitran_lines.temperatureDependence(w_index) .*...
    (hitran_lines.airBroadenedWidth(w_index) .* (P - P_self) +...
    hitran_lines.selfBroadenedWidth(w_index) .* P_self);


% compute the pressure broadened line shape for each line transition
f_pressure.shape = zeros(length(line_center), lineShape_length);
f_pressure.wavenum = zeros(length(line_center), lineShape_length);

for vv = 1:length(line_center)

    % Lines are narrow! Some less than a wavenumber
    f_pressure.wavenum(vv,:) = linspace(line_center(vv) - num_hwhm*f_pressure.hwhm(vv),...
        line_center(vv) + num_hwhm*f_pressure.hwhm(vv), lineShape_length);


    f_pressure.shape(vv,:) = f_pressure.hwhm(vv)./pi ./...
        (f_pressure.hwhm(vv).^2 + (f_pressure.wavenum(vv,:) - (line_center(vv) +...
        pressure_shift(vv) .* P)).^2 );

end



end