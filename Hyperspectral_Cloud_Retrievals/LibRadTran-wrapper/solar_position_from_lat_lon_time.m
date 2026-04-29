function [sza, phi0] = solar_position_from_lat_lon_time(date_in, time_utc_hours, lat, lon)
%SOLAR_POSITION_FROM_LAT_LON_TIME  Solar zenith and azimuth from date/time/lat/lon
%
%   [sza, phi0] = solar_position_from_lat_lon_time(date_in, time_utc_hours, lat, lon)
%
% Low-precision solar position algorithm based on the Astronomical Almanac
% (Michalsky 1988 / NOAA solar calculator). Accurate to ~0.01 deg between
% 1950-2050.
%
% INPUTS:
%   date_in        - MATLAB datetime (or string convertible to datetime).
%                    Year/month/day are taken from this; any time-of-day in
%                    date_in is ignored in favour of time_utc_hours.
%   time_utc_hours - scalar, decimal UTC time in hours (e.g., 14.5 for 14:30 UTC).
%   lat            - latitude in degrees (positive North).
%   lon            - longitude in degrees (positive East).
%
% OUTPUTS:
%   sza  - solar zenith angle in degrees (0 = sun directly overhead).
%   phi0 - solar azimuth angle in the libRadTran convention:
%          0-360 degrees clockwise from due south.

if ~isa(date_in, 'datetime')
    date_in = datetime(date_in);
end

% Build a UTC datetime that combines the calendar date with the decimal
% UTC hour-of-day value
dt = datetime(year(date_in), month(date_in), day(date_in)) + hours(time_utc_hours);

% Days since J2000.0 (2000 Jan 1.5 UT)
n = juliandate(dt) - 2451545.0;

% Mean longitude of the Sun (deg)
L = mod(280.460 + 0.9856474 * n, 360);
if L < 0, L = L + 360; end

% Mean anomaly (deg, then rad)
g = mod(357.528 + 0.9856003 * n, 360);
if g < 0, g = g + 360; end
g_rad = deg2rad(g);

% Ecliptic longitude (deg, then rad)
lambda     = L + 1.915 * sin(g_rad) + 0.020 * sin(2 * g_rad);
lambda_rad = deg2rad(lambda);

% Obliquity of the ecliptic (rad)
eps_rad = deg2rad(23.439 - 0.0000004 * n);

% Right ascension and declination (rad)
ra  = atan2(cos(eps_rad) * sin(lambda_rad), cos(lambda_rad));
dec = asin(sin(eps_rad) * sin(lambda_rad));

% Greenwich mean sidereal time (hours)
gmst = mod(6.697375 + 0.0657098242 * n + time_utc_hours, 24);
if gmst < 0, gmst = gmst + 24; end

% Local mean sidereal time (hours -> rad)
lmst = mod(gmst + lon / 15, 24);
if lmst < 0, lmst = lmst + 24; end
lmst_rad = deg2rad(lmst * 15);

% Hour angle (rad), wrapped to [-pi, pi]
ha = mod(lmst_rad - ra + pi, 2 * pi) - pi;

% Solar elevation and zenith (deg)
lat_rad = deg2rad(lat);
sin_el  = sin(dec) * sin(lat_rad) + cos(dec) * cos(lat_rad) * cos(ha);
el_rad  = asin(sin_el);
sza     = 90 - rad2deg(el_rad);

% Solar azimuth (rad) measured clockwise from due North
az_north_rad = atan2(-sin(ha), tan(dec) * cos(lat_rad) - sin(lat_rad) * cos(ha));
az_north     = mod(rad2deg(az_north_rad), 360);

% libRadTran convention: 0-360 deg clockwise from due South
phi0 = mod(az_north - 180, 360);

end
