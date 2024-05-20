%% Compute the Doppler broadened line shape (gaussian line shape)

% Doppler broadening tends to dominate in the upper atmosphere where
% pressure is low

% ----------- INPUTS -----------------

% (1) hitarn_lines - the hitran output as define by the function
% importhitran()

% (2) wavelength_grid - (nm) these wavelengths define the range over
% which the absorption cross section is desired. We must create line shapes for all
% absorption features within this range. *** MUST BE A CONTIGUOUS RANGE OF
% WAVELENGTHS, NOT DISPARATE VALUES ***

% (3) T - (K) Temperature of the gas


% ----------- OUTPUTS -----------------

% (1) f_doppler - ([cm^(-1)]^-1)a structure containing the independent variable, which
% is the wavenumber encompassing the doppler broadened line, the
% dependent variable, the Gaussian line shape, and the
% hwafl-width-at-half-max (HWHM)


% By Andrew John Buggee

%%

function f_doppler = gaussian_lineShape_for_hitran(hitran_lines, wavelength_grid, T)


%% Convert the wavelength grid to a wavenumber grid

% ------------------------------------------------------------------------
% ************ DOESN'T PLACE LINE SHAPES ON THE SAME GRID ****************
% ------------------------------------------------------------------------

% % Use the boundaries of this wavelength grid to make a finely spaced master
% % grid that all line shapes are computed on, such that their sum is made
% % easy
% d_lambda = 0.01;    % nm
% wavelength_master_grid = (wavelength_boundaries(1) - d_lambda*100):d_lambda:...
%     (wavelength_boundaries(end) + d_lambda*100);        % nm

% Make sure to convert the wavelength vector into microns
% Compute the line shapes on the master wavenumber grid so that we can
% easily sum over all line shapes to get a continuous curve of the
% absorption cross section.
wavenumber_master_grid = 10^4 ./ (wavelength_grid./1e3);        % cm^(-1)

% Make sure wavenumber_master_grid is a row vector
if size(wavenumber_master_grid,1)==1 && size(wavenumber_master_grid,2)>1

    wavenumber_master_grid = wavenumber_master_grid';
end

% Find the line centers closest to each value of the wavenumber grid
[~, w_index] = min(abs(wavenumber_master_grid - hitran_lines.transitionWavenumber'), [], 2);

% Drop any repeated indices
w_index = unique(w_index, 'stable');

% grab the line strength centers
line_center = hitran_lines.transitionWavenumber(w_index);


%% Compute the Half-Width-at-Half-Max and the Line Shape

% we have to compute a line shape for each spectral transition
lineShape_length = 100;

% define the number of HWHM that will define the total width of the
% wavenumber vector
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
f_doppler.hwhm = line_center./con.c .* sqrt((2*con.N_A * con.k_B * T * log(2))/molar_mass);      % cm^(-1)


% compute the Doppler line shape for each line transition

f_doppler.shape = zeros(length(line_center), lineShape_length);
f_doppler.wavenum = zeros(length(line_center), lineShape_length);

for vv = 1:length(line_center)

    % Lines are narrow! Some less than a wavenumber
    f_doppler.wavenum(vv,:) = linspace(line_center(vv) - num_hwhm*f_doppler.hwhm(vv),...
        line_center(vv) + num_hwhm*f_doppler.hwhm(vv), lineShape_length);

    f_doppler.shape(vv,:) = sqrt(log(2)./(pi*f_doppler.hwhm(vv).^2)) .* ...
        exp(- ((f_doppler.wavenum(vv,:) - line_center(vv)).^2 * log(2))./f_doppler.hwhm(vv)^2);

end



end