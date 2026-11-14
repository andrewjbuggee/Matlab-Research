%% Gaussian spectral response functions for the CPF HySICS channels (ver2)
%
% Same output as create_HySICS_specResponse.m when every channel has a
% 6-nm FWHM, without reading 'ILS Map.txt' and without fitting Gaussians to it.
%
% Why a ver2:
%   * Since 13 Oct 2026 the FWHM is fixed at 6 nm for every channel (Peter
%     Pilewskie, pers. comm.: the per-channel widths in 'ILS Map.txt' belong to
%     the balloon-tested HySICS instrument, not the CPF instrument). With a 
%     fixed FWHM only the center wavelengths are needed, so the ILS map read
%     (a readtable of the whole map, done on every call) and the 636 Gaussian
%     fits are dropped.
%   * The FWHM is an input (default 6 nm).
%   * Extra outputs for the one-run-per-cloud driver
%     (hysics_refl_from_synthetic_NN_inputs_ver2.m):
%       .center_wavelength  (n_bands x 1) nm, the Gaussian center of each band
%       .fwhm               (n_bands x 1) nm
%       .window             (n_bands x 2) nm, [first last] wavelength of the
%                           band's response grid = the range the band needs
%                           from a libRadtran run
%   * The center-wavelength file is closed after reading (the original leaves
%     the file open).
%
% The .value and .wavelength arrays are computed with exactly the same steps
% as create_HySICS_specResponse.m (center-wavelength extension, a ±20 nm grid
% on the source-file wavelength grid, a Gaussian evaluated on a fine grid and
% interpolated onto that grid). With fwhm_nm = 6 the two functions therefore
% return identical .value and .wavelength.
%
% INPUTS
%   band_number    vector of band numbers (1 ... 636)
%   sourceFile     libRadtran solar-source file name, e.g.
%                  'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat'
%                  (the response grid is defined on this file's wavelength grid;
%                  see create_HySICS_specResponse.m for the list of options)
%   computer_name  'anbu8374', 'andrewbuggee' or 'curc'
%   fwhm_nm        (optional) FWHM in nm for every channel. Default 6.
%
% OUTPUT
%   spec_response  struct with fields .value, .wavelength (n_bands x n_grid),
%                  .center_wavelength, .fwhm, .window, .fwhm_source
%
% By Andrew John Buggee
%%

function spec_response = create_HySICS_specResponse_ver2(band_number, sourceFile, computer_name, fwhm_nm)

if nargin < 3
    error([newline, 'Need at least three inputs: band numbers, source file, computer name', newline])
end
if nargin < 4 || isempty(fwhm_nm)
    fwhm_nm = 6;        % nm, CPF HySICS (P. Pilewskie, 13 Oct 2026)
end

if ~isnumeric(band_number)
    error('band_number must be numeric')
end
if any(band_number < 1) || any(band_number > 636)
    error('HySICS has 636 spectral bands; band_number must lie in 1 ... 636')
end


% --- folder holding 'Center Wavelength.txt' ---------------------------------
if strcmp(computer_name, 'anbu8374')
    folderPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/';
elseif strcmp(computer_name, 'andrewbuggee')
    folderPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/';
elseif strcmp(computer_name, 'curc')
    folderPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/';
else
    error('Unknown computer name: %s', computer_name)
end


% --- center wavelengths of the channels in the instrument file (407-1813 nm) ---
file_id = fopen([folderPath, 'Center Wavelength.txt'], 'r');
if file_id < 0
    error('Could not open %s', [folderPath, 'Center Wavelength.txt'])
end
wl_center = textscan(file_id, '%f', 'Delimiter', ',', 'MultipleDelimsAsOne', 1, ...
    'CommentStyle', ';');
fclose(file_id);
wl_center = wl_center{1};

% Extend to HySICS' full 350-2300 nm range with the mean channel spacing
% (identical to create_HySICS_specResponse.m)
mean_wl_diff = mean(diff(wl_center));        % nm
full_center_wl = [fliplr(wl_center(1)-mean_wl_diff:-mean_wl_diff:350)'; wl_center; ...
    (wl_center(end)+mean_wl_diff:mean_wl_diff:2300)'];


% --- source-file wavelength grid -----------------------------------------------
% libRadtran reports radiance on the wavelength grid of the solar source file,
% so each response function is defined on that grid.
[~, source_wavelength] = read_solar_flux_file([full_center_wl(1)-25, full_center_wl(end)+25], sourceFile);  % nm
sourceFile_wavelength_resolution = source_wavelength(2) - source_wavelength(1);   % nm
source_wavelength_grid = -20:sourceFile_wavelength_resolution:20;                 % nm, ±20 nm window

% fine grid on which the continuous Gaussian is evaluated before interpolation
specResponse_wavelength_grid = -21:(sourceFile_wavelength_resolution/20):21;      % nm

fwhm_all = fwhm_nm * ones(length(full_center_wl), 1);    % nm, same for every channel


% --- Gaussian response of each requested band ----------------------------------
n_bands = length(band_number);
spec_response.value      = zeros(n_bands, length(source_wavelength_grid));
spec_response.wavelength = zeros(n_bands, length(source_wavelength_grid));
spec_response.center_wavelength = zeros(n_bands, 1);
spec_response.fwhm              = zeros(n_bands, 1);

for nn = 1:n_bands

    bb = band_number(nn);

    % standard deviation of a Gaussian from its FWHM: FWHM = 2*sqrt(2*ln 2)*sigma
    sigma = fwhm_all(bb) / (2*sqrt(2*log(2)));

    % high-resolution grid centerd on the channel's center wavelength
    wl_grid_hiRes = specResponse_wavelength_grid + full_center_wl(bb);

    % the output grid is the ±20 nm source-grid window centerd on the source
    % wavelength closest to the channel center (so it is NOT exactly centerd on
    % the Gaussian peak)
    [~, center_wavelength_idx] = min(abs(source_wavelength - full_center_wl(bb)));
    spec_response.wavelength(nn, :) = source_wavelength_grid + source_wavelength(center_wavelength_idx);

    % Gaussian on the fine grid, interpolated to the source grid
    response_hiRes = pdf('Normal', wl_grid_hiRes, full_center_wl(bb), sigma);
    spec_response.value(nn, :) = interp1(wl_grid_hiRes, response_hiRes, spec_response.wavelength(nn, :));

    spec_response.center_wavelength(nn) = full_center_wl(bb);
    spec_response.fwhm(nn)              = fwhm_all(bb);

end

spec_response.window = [spec_response.wavelength(:, 1), spec_response.wavelength(:, end)];   % nm
spec_response.fwhm_source = sprintf('fixed FWHM = %g nm for every channel', fwhm_nm);

end
