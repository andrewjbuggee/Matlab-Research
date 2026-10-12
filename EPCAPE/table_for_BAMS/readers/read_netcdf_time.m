function t = read_netcdf_time(file, var_name)
% READ_NETCDF_TIME  Read a netCDF time variable as a MATLAB datetime (UTC).
%
%   t = read_netcdf_time(file)              reads the variable 'time'
%   t = read_netcdf_time(file, var_name)
%
%   netCDF files store time as a number with a units attribute such as
%     "seconds since 1970-01-01"                (the combined ARM files)
%     "seconds since 2023-02-15"                (each SONDEPARAM launch file)
%     "nanoseconds since 2023-02-15 00:00:00"   (Abbey Williams' 15-min file)
%   This function reads the units, finds the reference date after "since",
%   and adds the stored offsets to it.
%
%   Output   t   datetime column vector (no time zone attached; all times
%                in these files are UTC)

if nargin < 2
    var_name = 'time';
end

offsets = double(ncread(file, var_name));          % stored numbers
units   = strtrim(string(ncreadatt(file, var_name, 'units')));

% split "<unit> since <reference date>"
parts = split(units, " since ");
unit_name = lower(strtrim(parts(1)));
ref_text  = strtrim(parts(2));
ref_text  = erase(ref_text, " UTC");               % e.g. "1970-01-01 00:00:00 UTC"
ref_text  = erase(ref_text, " 0:00");              % e.g. "2023-02-15 00:00:00 0:00"
if strlength(ref_text) == 10                       % date only, e.g. "2023-02-15"
    ref_time = datetime(ref_text, 'InputFormat', 'yyyy-MM-dd');
else
    ref_time = datetime(ref_text, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
end

switch unit_name
    case "seconds"
        t = ref_time + seconds(offsets);
    case "nanoseconds"
        % round to whole seconds: the 15-min file has a few nanoseconds of jitter
        t = ref_time + seconds(round(offsets / 1e9));
    otherwise
        error('read_netcdf_time: unexpected time units "%s" in %s', units, file);
end
t = t(:);

end
