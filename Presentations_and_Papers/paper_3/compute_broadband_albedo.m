function albedo = compute_broadband_albedo(re, lwc, z_km, RT, folder_paths, ...
    index, profile_tag, date_str, time_val, campaign_name)
%% COMPUTE_BROADBAND_ALBEDO
%
% Writes a water-cloud file and libRadtran input file, runs uvspec, reads
% the integrated solar-flux output, and returns the broadband (250–4000 nm)
% shortwave cloud albedo at the top of atmosphere.
%
% Albedo = E_up(TOA) / E_dir(TOA)
%   where E_dir(TOA) = F0 * cos(SZA) — the direct solar beam at TOA
%   and   E_up(TOA)  = the diffuse upwelling flux at TOA
%
% INPUTS
%   re           – effective radius profile [µm]  (column vector, top→bottom)
%   lwc          – liquid water content profile [g/m³] (column vector, top→bottom)
%   z_km         – altitude profile [km MSL]  (column vector, top→bottom)
%   RT           – structure with libRadtran settings:
%                    .sza            solar zenith angle [degrees]
%                    .phi0           solar azimuth angle [degrees]
%                    .day_of_year    integer day of year
%                    .surface_albedo scalar surface albedo
%                    .atm_file       atmosphere file name (e.g. 'afglus.dat')
%                    .source_file    solar source file name
%   folder_paths – structure from define_folderPaths_for_HySICS:
%                    .libRadtran_inp
%                    .libRadtran_data
%                    .libRadtran_water_cloud_files
%                    .which_computer
%   index        – integer, unique identifier for temporary files
%   profile_tag  – short string label used in filename (e.g. 'insitu')
%   date_str     – date string for wc-file naming
%   time_val     – UTC time for wc-file naming
%   campaign_name – 'vocalsRex' or 'oracles'
%
% OUTPUT
%   albedo  – broadband shortwave albedo [0–1]
%
% By Andrew John Buggee

%% Input validation

if nargin ~= 10
    error([mfilename, ': requires 10 inputs.'])
end

if length(re) ~= length(lwc) || length(re) ~= length(z_km)
    error([mfilename, ': re, lwc, z_km must be the same length.'])
end

if any(isnan(re)) || any(re <= 0)
    error([mfilename, ': re contains NaN or non-positive values.'])
end

if any(lwc < 0)
    error([mfilename, ': lwc contains negative values.'])
end

%% Write the water-cloud file

% Use a unique filename combining profile_tag and index
unique_index = index;   % re-use the calling-loop index

wc_filename = write_wc_file_from_in_situ(re(:), lwc(:), z_km(:), ...
    campaign_name, date_str, time_val, ...
    false, ...                              % compute_weighting_functions = false
    folder_paths.which_computer, ...
    unique_index, ...
    folder_paths.libRadtran_water_cloud_files);

%% Write the libRadtran INP file for broadband solar flux

inp_basename = sprintf('broadband_albedo_%s_%04d', profile_tag, unique_index);
inp_name     = [inp_basename, '.INP'];
out_name     = [inp_basename, '.OUT'];

write_broadband_INP(folder_paths.libRadtran_inp, ...
    folder_paths.libRadtran_data, ...
    folder_paths.libRadtran_water_cloud_files, ...
    inp_name, wc_filename, RT);

%% Run uvspec

runUVSPEC_ver2(folder_paths.libRadtran_inp, inp_name, inp_basename, ...
    folder_paths.which_computer);

%% Read the output and compute albedo

out_file = [folder_paths.libRadtran_inp, out_name];

if ~exist(out_file, 'file')
    error([mfilename, ': uvspec output file not found:\n  %s'], out_file)
end

% Expected output format (one line, output_process integrate):
%   edir  edn  eup  uavgdir  uavgdn  uavgup
out_data = load(out_file);

if isempty(out_data) || size(out_data, 2) < 3
    error([mfilename, ': uvspec output file empty or missing columns:\n  %s'], out_file)
end

% Column layout (see libRadtran manual, §4.5):
%   1: edir   — direct beam irradiance at zout level [W/m²]
%   2: edn    — diffuse downwelling irradiance [W/m²]
%   3: eup    — diffuse upwelling irradiance [W/m²]
% (when output_process integrate is used, wavelength column is dropped)

edir = out_data(1);   % W/m² integrated
edn  = out_data(2);   % W/m²  (= 0 at TOA — nothing above)
eup  = out_data(3);   % W/m²  (reflected solar radiation)

% Broadband albedo = reflected / incoming
incoming = edir + edn;   % Total downwelling at TOA
if incoming < 1
    % Protect against degenerate cases (e.g. SZA = 90°)
    albedo = NaN;
    warning([mfilename, ': total incoming flux = %.3f W/m² (profile %d, %s).'], ...
        incoming, index, profile_tag)
    return
end

albedo = eup / incoming;

% Sanity check
if albedo < 0 || albedo > 1
    warning([mfilename, ': albedo = %.4f out of [0,1] range (profile %d, %s).'], ...
        albedo, index, profile_tag)
    albedo = max(0, min(1, albedo));
end

end
