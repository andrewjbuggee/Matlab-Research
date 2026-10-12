function F = read_15min_file(path)
% READ_15MIN_FILE  Read the variables used here from Abbey Williams' 15-min file.
%
%   F = read_15min_file(path)          path to EPCAPE_15mindata.nc
%
%   EPCAPE_15mindata.nc (Abbey Williams) merges many EPCAPE records onto one
%   15-min time axis, 2023-02-15 00:00 to 2024-02-15 00:00 (UTC), each
%   interval labelled by its START. This is the file Kavin's AMSforAndrew.mlx
%   uses for the cloud/haze rows.
%
%   Output: struct of column vectors (one element per 15-min interval)
%     F.time        interval start (datetime, UTC)
%     F.class       ClrCldHaz_class: 1 = cloud, 2 = haze, 0 = CVI off,
%                   NaN = CVI changing state. Cloud and haze are only set
%                   while the CVI inlet was on (CVIon_class = 1).
%     F.nd_cm3      FM120_Nd  : FM-120 droplet number, 2-50 um (# cm-3)
%     F.lwc_gm3     FM120_LWC : FM-120 liquid water content (g m-3)
%     F.dsd         FM120_DSD : droplets per size bin (# cm-3), matrix (intervals x 30 bins)
%     F.diam_um     FM120_diam: the diameter given to each bin (um). These are
%                   the bins' UPPER edges (3, 4, ..., 50 um for the 2-3, ...,
%                   48-50 um bins listed in readme_FM120.md).
%     F.ams_ugm3    AMS_mass  : AMS mass concentration (ug m-3), matrix (intervals x 5)
%     F.species     names of the 5 AMS columns, in file order:
%                   nitrate, sulfate, organics, ammonium, chloride
%                   While the CVI is on, these are the CE-corrected AMS values
%                   divided by the GCVI enhancement factor (droplet residuals).

F.time  = read_netcdf_time(path, 'time');
F.class = double(ncread(path, 'ClrCldHaz_class'));
F.nd_cm3  = double(ncread(path, 'FM120_Nd'));
F.lwc_gm3 = double(ncread(path, 'FM120_LWC'));

% netCDF stores these as (time, bin); MATLAB's ncread returns the dimensions
% in the reverse order, (bin, time), so transpose to (time, bin).
F.dsd      = double(ncread(path, 'FM120_DSD')).';
F.diam_um  = double(ncread(path, 'FM120_diam'));
F.ams_ugm3 = double(ncread(path, 'AMS_mass')).';

F.species = lower(strtrim(string(ncread(path, 'AMS_species'))));
expected  = ["nitrate"; "sulfate"; "organics"; "ammonium"; "chloride"];
if ~isequal(F.species(:), expected)
    error('read_15min_file: unexpected AMS_species order: %s', strjoin(F.species, ', '));
end

end
