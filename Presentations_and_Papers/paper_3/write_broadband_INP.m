function write_broadband_INP(inp_folder, data_folder, wc_folder, ...
    inp_filename, wc_filename, RT)
%% WRITE_BROADBAND_INP
%
% Writes a minimal libRadtran .INP file for computing broadband (250–4000 nm)
% solar irradiance fluxes at the top of atmosphere.
%
% Output is the integrated (250–4000 nm) solar flux using
%   output_process integrate
%
% Output columns (edir  edn  eup) at TOA:
%   edir — direct beam downwelling irradiance [W/m²]
%   edn  — diffuse downwelling irradiance [W/m²]
%   eup  — diffuse upwelling irradiance [W/m²]
%
% Broadband albedo = eup / (edir + edn)
%
% INPUTS
%   inp_folder    – folder where .INP file will be written
%   data_folder   – libRadtran data folder (ends with '/')
%   wc_folder     – water-cloud files folder (ends with '/')
%   inp_filename  – output .INP filename (e.g. 'broadband_albedo_insitu_0001.INP')
%   wc_filename   – water-cloud .DAT filename (just the filename, not full path)
%   RT            – structure:
%                     .sza            solar zenith angle [degrees]
%                     .phi0           solar azimuth angle [degrees]
%                     .day_of_year    integer day of year
%                     .surface_albedo scalar
%                     .atm_file       atmosphere file (e.g. 'afglus.dat')
%                     .source_file    solar source file
%
% By Andrew John Buggee

%% Open file

fid = fopen([inp_folder, inp_filename], 'w');
if fid < 0
    error('write_broadband_INP: could not open file for writing:\n  %s', ...
        [inp_folder, inp_filename])
end

%% RTE solver

fprintf(fid, 'rte_solver disort\n\n');

%% Number of streams

fprintf(fid, 'number_of_streams 16\n\n');

%% Nakajima–Tanaka radiance correction

fprintf(fid, 'disort_intcor phase    # Nakajima & Tanaka phase function correction\n\n');

%% Data files path

fprintf(fid, 'data_files_path %s\n\n', data_folder);

%% Band parameterisation  (reptran coarse — balanced speed / accuracy)

fprintf(fid, 'mol_abs_param reptran coarse\n\n');

%% Atmosphere file

fprintf(fid, 'atmosphere_file %satmmod/%s\n\n', data_folder, RT.atm_file);

%% Solar source

fprintf(fid, 'source solar %s\n\n', RT.source_file);

%% Day of year (Earth–Sun distance correction)

fprintf(fid, 'day_of_year %d\n\n', RT.day_of_year);

%% Surface albedo

fprintf(fid, 'albedo %.4f    # Ocean surface albedo\n\n', RT.surface_albedo);

%% Water cloud file

fprintf(fid, 'wc_file 1D %s%s\n', wc_folder, wc_filename);
fprintf(fid, 'wc_properties mie interpolate\n\n');

%% Solar zenith and azimuth angles

fprintf(fid, 'sza %.2f    # Solar zenith angle [degrees]\n\n', RT.sza);
fprintf(fid, 'phi0 %.2f    # Solar azimuth angle [degrees from due South]\n\n', RT.phi0);

%% Wavelength range — full solar spectrum

fprintf(fid, 'wavelength 250 4000    # nm: full solar spectral range\n\n');

%% Output level: top of atmosphere only

fprintf(fid, 'zout toa    # Output at top of atmosphere\n\n');

%% Output quantities: direct, diffuse-down, diffuse-up fluxes

fprintf(fid, 'output_user edir edn eup    # W/m^2/nm at each wavelength\n\n');

%% Integrate over wavelength

fprintf(fid, 'output_process integrate    # Produce one integrated value (W/m^2)\n\n');

%% Quiet error messages

fprintf(fid, 'quiet\n');

%% Close

fclose(fid);

end
