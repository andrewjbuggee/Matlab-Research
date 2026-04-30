%% Write a libRadtran .DAT radiosonde file from a synthetic ERA5-style profile
%
% This is the synthetic-cloud counterpart to write_ERA5_radiosonde_DAT_with_multiPixels.
% Instead of looking up the ERA5 grid cell nearest a real (lat, lon, datetime),
% it takes pressure, temperature, and water-vapor number density already in
% memory (typically read from the per-cloud row of the training-inputs .nc file
% produced by 07_build_training_inputs.py) and writes them out in libRadtran
% radiosonde format.
%
% INPUTS:
%   pressure_hPa       (n_levels x 1) Pressure on the ERA5 grid                [hPa]
%   temperature_K      (n_levels x 1) Temperature                              [K]
%   vapor_concentration (n_levels x 1) Water vapor number density              [molec/cm^3]
%   folder_paths       struct with field .atm_folder_path (where to write)
%   filename           (optional) custom file name. Pass [] to auto-generate
%                      using the cloud index. Must end with .dat or .DAT.
%   num_vars           2 -> write only p,T;  3 -> write p,T,vapor (default 3)
%   cloud_idx          (optional) integer index used in the auto-generated
%                      filename. Pass [] to omit.
%   print_status_updates  (optional, default true)
%
% OUTPUTS:
%   filename_fullPath  Full path to the written .DAT file. Pass this string
%                      back as inputs.RT.radiosonde_file in the libRadtran INP.
%
% Notes:
%   - libRadtran wants the radiosonde rows ordered TOA -> surface
%     (pressure increasing).  This routine reorders if needed.
%   - No NaN handling beyond a warning; synthetic profiles should be clean.
%   - No specific-humidity → number-density conversion: the synthetic
%     pipeline already stores water vapor as molec/cm^3.
%
% By Andrew John Buggee  (synthetic-input variant, drafted with Claude)
%%

function filename_fullPath = write_synthetic_radiosonde_DAT(pressure_hPa, ...
        temperature_K, vapor_concentration, folder_paths, filename, ...
        num_vars, cloud_idx, print_status_updates)

% --- Defaults ----------------------------------------------------------------
if nargin < 6 || isempty(num_vars)
    num_vars = 3;
end
if nargin < 7
    cloud_idx = [];
end
if nargin < 8 || isempty(print_status_updates)
    print_status_updates = false;
end

if ~isfield(folder_paths, 'atm_folder_path')
    error('folder_paths struct must contain field "atm_folder_path".');
end

atm_folder_path = folder_paths.atm_folder_path;
if atm_folder_path(end) ~= filesep && atm_folder_path(end) ~= '/'
    atm_folder_path = [atm_folder_path, '/'];
end

% --- Coerce to column vectors -----------------------------------------------
pressure_hPa        = double(pressure_hPa(:));
temperature_K       = double(temperature_K(:));
vapor_concentration = double(vapor_concentration(:));

if ~(numel(pressure_hPa) == numel(temperature_K) ...
     && numel(temperature_K) == numel(vapor_concentration))
    error('pressure / temperature / vapor must be the same length (got %d, %d, %d).', ...
        numel(pressure_hPa), numel(temperature_K), numel(vapor_concentration));
end

if any(isnan(temperature_K)) || any(isnan(vapor_concentration)) || any(isnan(pressure_hPa))
    warning('Synthetic radiosonde profile contains NaN values.');
end

% --- libRadtran wants TOA -> surface (pressure increasing) -------------------
if pressure_hPa(1) > pressure_hPa(end)
    p_w   = flipud(pressure_hPa);
    T_w   = flipud(temperature_K);
    vap_w = flipud(vapor_concentration);
else
    p_w   = pressure_hPa;
    T_w   = temperature_K;
    vap_w = vapor_concentration;
end

% --- Filename ---------------------------------------------------------------
if isempty(filename)
    if isempty(cloud_idx)
        filename_radiosonde = 'synthetic_radiosonde.dat';
    else
        filename_radiosonde = sprintf('synthetic_radiosonde_cloud%05d.dat', cloud_idx);
    end
else
    filename_radiosonde = filename;
    if ~contains(filename_radiosonde, '.dat') && ~contains(filename_radiosonde, '.DAT')
        filename_radiosonde = [filename_radiosonde, '.dat'];
    end
end

filename_fullPath = [atm_folder_path, filename_radiosonde];

% --- Write the file ---------------------------------------------------------
fid = fopen(filename_fullPath, 'w');
if fid == -1
    error('Could not open file for writing: %s', filename_fullPath);
end

if num_vars == 2
    fprintf(fid, '# synthetic ERA5-style profile (T, P only)\n');
    fprintf(fid, '#   p(hPa)  T(K)\n');
    for ii = 1:length(p_w)
        fprintf(fid, '%12.5f %5.1f\n', p_w(ii), T_w(ii));
    end
elseif num_vars == 3
    fprintf(fid, '# synthetic ERA5-style profile (T, P, water-vapor number density)\n');
    fprintf(fid, '#   p(hPa)  T(K)  H2O(cm-3)\n');
    for ii = 1:length(p_w)
        fprintf(fid, '%f %f %e\n', p_w(ii), T_w(ii), vap_w(ii));
    end
else
    fclose(fid);
    error('num_vars must be 2 or 3, got %d', num_vars);
end

fclose(fid);

if print_status_updates
    fprintf('synthetic radiosonde written: %s\n', filename_fullPath);
    fprintf('  p   range: %.2f .. %.2f hPa\n', min(p_w), max(p_w));
    fprintf('  T   range: %.2f .. %.2f K\n',   min(T_w), max(T_w));
    if num_vars == 3
        fprintf('  vap range: %.3e .. %.3e molec/cm^3\n', min(vap_w), max(vap_w));
    end
end

end
