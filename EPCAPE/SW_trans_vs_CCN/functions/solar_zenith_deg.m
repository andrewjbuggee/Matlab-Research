%% Solar zenith angle at a site, from the NOAA general solar position equations
%
% DESCRIPTION:
%   Solar zenith angle (deg) from the NOAA Global Monitoring Laboratory
%   "General Solar Position Calculations": the fractional year gives the
%   equation of time and the solar declination (Fourier series of Spencer
%   1971), then the hour angle from the true solar time. Accuracy ~0.1-0.2
%   deg, enough for day/night shading (refraction is not included; use
%   SZA > 90.833 deg for "sun below the horizon" as NOAA does for sunrise
%   and sunset).
%
% INPUTS:
%   t_utc   - times (datetime, UTC)
%   lat_deg - latitude (deg north)
%   lon_deg - longitude (deg east; negative = west)
%
% OUTPUTS:
%   sza_deg - solar zenith angle (deg), same size as t_utc
%
% REFERENCES:
%   NOAA Global Monitoring Laboratory, "General Solar Position
%     Calculations", https://gml.noaa.gov/grad/solcalc/solareqns.PDF.
%   Spencer, J. W., 1971: Fourier series representation of the position of
%     the Sun. Search, 2, 172.
%
% By Andrew John Buggee
%%

function sza_deg = solar_zenith_deg(t_utc, lat_deg, lon_deg)

doy  = day(t_utc, 'dayofyear');
hr   = hour(t_utc) + minute(t_utc) / 60 + second(t_utc) / 3600;   % UTC hour of day
n_yr = 365 + (mod(year(t_utc), 4) == 0);                          % days in the year (valid to 2099)
g = 2 * pi ./ n_yr .* (doy - 1 + (hr - 12) / 24);                 % fractional year (rad)

% equation of time (min) and solar declination (rad)
eqtime = 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) ...
                   - 0.014615 * cos(2 * g) - 0.040849 * sin(2 * g));
decl = 0.006918 - 0.399912 * cos(g) + 0.070257 * sin(g) - 0.006758 * cos(2 * g) ...
       + 0.000907 * sin(2 * g) - 0.002697 * cos(3 * g) + 0.00148 * sin(3 * g);

% true solar time (min; time zone 0 because t is UTC) and hour angle (deg)
tst = hr * 60 + eqtime + 4 * lon_deg;
ha_deg = tst / 4 - 180;

cos_z = sind(lat_deg) .* sin(decl) + cosd(lat_deg) .* cos(decl) .* cosd(ha_deg);
sza_deg = acosd(min(max(cos_z, -1), 1));

end
