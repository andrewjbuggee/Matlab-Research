%% Read ORACLES Microphysics Dataset (combined_microphysics netCDF files)

% INPUTS:
% -------
%   filename - full path to the ORACLES combined microphysics netCDF file
%              (e.g., 'Microphysics_P3_YYYYMMDD_R0.nc')
%
% OUTPUTS:
% --------
%   oracles  - structure with fields analogous to readVocalsRex output:
%
%     Instruments:
%       CAS    : 3  < D < 45 um (analogous to CDP in VOCALS-REx)
%       2DS    : 45 < D < 975 um (analogous to 2DC)
%       HVPS-3 : 975 < D < 10575 um
%
%     Size distribution:
%       Nd                   - (N_bins x N_time) number concentration per bin [cm^-3]
%       drop_radius_bin_center - (1 x 172) radius bin centers [µm]
%       drop_radius_bin_edges  - (1 x 173) radius bin edges [µm]
%
%     Total quantities (all probes):
%       total_Nc    - total droplet number concentration [cm^-3]
%       re          - effective radius (pre-computed, all probes) [µm]
%       lwc         - liquid water content (all probes) [g/m^3]
%
%     CAS probe (D <= 50 µm, analogous to CDP):
%       total_Nc_CAS - cloud droplet number concentration [cm^-3]
%       re_CAS       - effective radius from CAS bins [µm]
%       cwc          - cloud water content (D <= 50 µm) [g/m^3]
%
%     2DS+HVPS-3 probes (D > 50 µm, analogous to 2DC):
%       total_Nc_2DS_HVPS - drizzle/rain droplet concentration [cm^-3]
%       re_2DS_HVPS       - effective radius from 2DS+HVPS bins [µm]
%       rwc               - rain water content (D > 50 µm) [g/m^3]
%
%     Other cloud properties:
%       bulk_lwc    - King hot-wire LWC [g/m^3]
%       r_rain      - rain rate [mm/hr]
%       Na          - out-of-cloud aerosol conc (PCASP) [cm^-3]
%       Z           - radar reflectivity [dBZ]
%       beta        - extinction [km^-1]
%       rv          - mean volume radius [µm]
%
%     Navigation and atmosphere:
%       time        - seconds since midnight UTC [s]
%       timevec     - MATLAB datenum [days]
%       dateOfFlight - datetime of the flight
%       latitude    - [degrees North]
%       longitude   - [degrees East]
%       altitude    - [m above MSL]
%       temp        - static air temperature [degrees C]
%       pres        - static pressure [millibars]

% By Andrew John Buggee

%%

function oracles = readORACLES(filename)


% -------------------------------------------------------------------------
% ---------------------- Flight date and time -----------------------------
% -------------------------------------------------------------------------

% Extract date from filename: 'Microphysics_P3_YYYYMMDD_R0.nc'
[~, fname, ~] = fileparts(filename);
date_str = fname(17:24);        % Extract YYYYMMDD
dateOfFlight = datetime(date_str, 'InputFormat', 'yyyyMMdd');

% Read the timevec variable (MATLAB datenum: days since Jan 0, 0000)
timevec = double(ncread(filename, 'timevec'));           % days
timevec = reshape(timevec, 1, []);

% Convert to seconds since midnight UTC
time = (mod(timevec, 1)) * 86400;                       % seconds since midnight UTC
time = reshape(time, 1, []);

N_time = length(time);


% -------------------------------------------------------------------------
% ----------------------- Size bin information ----------------------------
% -------------------------------------------------------------------------
% Bins are defined in diameter space (µm); convert to radius

bin_min = double(ncread(filename, 'bin_min'));           % µm diameter, lower edge (172x1)
bin_max = double(ncread(filename, 'bin_max'));           % µm diameter, upper edge (172x1)
bin_mid = double(ncread(filename, 'bin_mid'));           % µm diameter, center (172x1)

% Bins are contiguous: bin_max(i) == bin_min(i+1)
% Construct edge array from min values plus the final max
drop_diameter_bin_edges = [bin_min; bin_max(end)];          % 173x1, µm diameter
drop_radius_bin_edges   = reshape(drop_diameter_bin_edges / 2, 1, []);  % 1x173, µm radius
drop_radius_bin_center  = reshape(bin_mid / 2, 1, []);                  % 1x172, µm radius

% Index for CAS-equivalent bins (D <= 50 µm, analogous to CDP)
% This matches the ORACLES Nc variable definition (D < 50 µm)
idx_CAS = (bin_mid <= 50)';     % 172x1 logical


% -------------------------------------------------------------------------
% -------------------- Size distribution (Nd) -----------------------------
% -------------------------------------------------------------------------
% Nd is (N_bins x N_time) = (172 x N_time), units: cm^-3 per bin
% Outside-cloud values are NaN; replace with 0 for threshold comparisons

Nd = double(ncread(filename, 'Nd'))';                    % cm^-3 per bin (172 x N_time)
Nd(isnan(Nd)) = 0;


% -------------------------------------------------------------------------
% ----------- Compute effective radii from size distribution --------------
% -------------------------------------------------------------------------
% r^3 / r^2 moment ratio, analogous to vocalsRex approach

% Radius matrix in cm (for unit consistency)
r_matrix_cm = repmat(drop_radius_bin_center', 1, N_time) ./ 1e4;   % cm (172 x N_time)

% Effective radius from CAS bins (D <= 50 µm)
re_CAS = double(sum(r_matrix_cm(idx_CAS,:).^3 .* Nd(idx_CAS,:), 1) ./ ...
    sum(r_matrix_cm(idx_CAS,:).^2 .* Nd(idx_CAS,:), 1) * 1e4);    % µm
re_CAS(isnan(re_CAS)) = 0;

% Effective radius from 2DS+HVPS-3 bins (D > 50 µm)
idx_2DS_HVPS = ~idx_CAS;
re_2DS_HVPS = double(sum(r_matrix_cm(idx_2DS_HVPS,:).^3 .* Nd(idx_2DS_HVPS,:), 1) ./ ...
    sum(r_matrix_cm(idx_2DS_HVPS,:).^2 .* Nd(idx_2DS_HVPS,:), 1) * 1e4);  % µm
re_2DS_HVPS(isnan(re_2DS_HVPS)) = 0;


% -------------------------------------------------------------------------
% ----------- Pre-computed cloud microphysical properties -----------------
% -------------------------------------------------------------------------
% These are quality-controlled values from the ORACLES data providers.
% Outside-cloud values are NaN; replace with 0.

% Total (all probes)
total_Nc    = reshape(double(ncread(filename, 'N')),  1, []);   % cm^-3
Re          = reshape(double(ncread(filename, 'Re')), 1, []);   % µm
lwc         = reshape(double(ncread(filename, 'lwc')),1, []);   % g/m^3

% CAS-equivalent (D <= 50 µm)
total_Nc_CAS = reshape(double(ncread(filename, 'Nc')),  1, []); % cm^-3
cwc          = reshape(double(ncread(filename, 'cwc')), 1, []); % g/m^3

% 2DS+HVPS-equivalent (D > 50 µm)
total_Nc_2DS_HVPS = reshape(double(ncread(filename, 'Nc50')), 1, []); % cm^-3
rwc               = reshape(double(ncread(filename, 'rwc')),  1, []); % g/m^3

% Other cloud properties
bulk_lwc = reshape(double(ncread(filename, 'bulk_lwc')), 1, []); % g/m^3 (King hot-wire)
r_rain   = reshape(double(ncread(filename, 'r')),        1, []); % mm/hr
Na       = reshape(double(ncread(filename, 'Na')),       1, []); % cm^-3
Z        = reshape(double(ncread(filename, 'Z')),        1, []); % dBZ
beta     = reshape(double(ncread(filename, 'beta')),     1, []); % km^-1
rv       = reshape(double(ncread(filename, 'rv')),       1, []); % µm

% Replace NaN (and Inf for lwc) with 0 for out-of-cloud regions
total_Nc(isnan(total_Nc))                     = 0;
Re(isnan(Re))                                 = 0;
lwc(isnan(lwc) | isinf(lwc))                  = 0;
total_Nc_CAS(isnan(total_Nc_CAS))             = 0;
cwc(isnan(cwc))                               = 0;
total_Nc_2DS_HVPS(isnan(total_Nc_2DS_HVPS))   = 0;
rwc(isnan(rwc))                               = 0;
bulk_lwc(isnan(bulk_lwc))                     = 0;
r_rain(isnan(r_rain))                         = 0;
Na(isnan(Na))                                 = 0;
Z(isnan(Z))                                   = 0;
beta(isnan(beta))                             = 0;
rv(isnan(rv))                                 = 0;


% -------------------------------------------------------------------------
% ----------- Navigation and atmospheric state variables ------------------
% -------------------------------------------------------------------------

altitude  = reshape(double(ncread(filename, 'altitude')),  1, []); % m MSL
latitude  = reshape(double(ncread(filename, 'latitude')),  1, []); % degrees N
longitude = reshape(double(ncread(filename, 'longitude')), 1, []); % degrees E
temp      = reshape(double(ncread(filename, 'temp')),      1, []); % degrees C
pres      = reshape(double(ncread(filename, 'pres')),      1, []); % millibars


% -------------------------------------------------------------------------
% ----------------------- Assemble output structure -----------------------
% -------------------------------------------------------------------------

% Size distribution and bin information
oracles.Nd                     = Nd;                        % cm^-3 per bin (172 x N_time)
oracles.drop_radius_bin_center = drop_radius_bin_center;    % µm radius, 1x172
oracles.drop_radius_bin_edges  = drop_radius_bin_edges;     % µm radius, 1x173

% Total number concentration and LWC (all probes combined)
oracles.total_Nc    = total_Nc;     % cm^-3
oracles.re          = Re;           % µm (pre-computed by ORACLES data group)
oracles.lwc         = lwc;          % g/m^3

% CAS probe quantities (D <= 50 µm, analogous to CDP in VOCALS-REx)
oracles.total_Nc_CAS = total_Nc_CAS;    % cm^-3
oracles.re_CAS       = re_CAS;          % µm (computed from CAS bins)
oracles.cwc          = cwc;             % g/m^3 (cloud water content)

% 2DS+HVPS-3 probe quantities (D > 50 µm, analogous to 2DC in VOCALS-REx)
oracles.total_Nc_2DS_HVPS = total_Nc_2DS_HVPS; % cm^-3
oracles.re_2DS_HVPS       = re_2DS_HVPS;        % µm (computed from 2DS+HVPS bins)
oracles.rwc               = rwc;                % g/m^3 (rain water content)

% Other cloud microphysical properties
oracles.bulk_lwc  = bulk_lwc;   % g/m^3 (King hot-wire LWC)
oracles.r_rain    = r_rain;     % mm/hr (rain rate)
oracles.Na        = Na;         % cm^-3 (out-of-cloud PCASP aerosol)
oracles.Z         = Z;          % dBZ (radar reflectivity)
oracles.beta      = beta;       % km^-1 (extinction)
oracles.rv        = rv;         % µm (mean volume radius)

% Time and navigation
oracles.time         = time;            % s since midnight UTC
oracles.timevec      = timevec;         % MATLAB datenum [days]
oracles.dateOfFlight = dateOfFlight;    % datetime
oracles.latitude     = latitude;        % degrees N
oracles.longitude    = longitude;       % degrees E
oracles.altitude     = altitude;        % m MSL

% Atmospheric state
oracles.temp = temp;    % degrees C
oracles.pres = pres;    % millibars


end
