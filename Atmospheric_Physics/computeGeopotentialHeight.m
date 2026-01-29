% computeGeopotentialHeight - Calculate geopotential height from ERA5 data
%
% DESCRIPTION:
%   Computes geopotential height at pressure levels using the hypsometric
%   equation (also called the thickness equation). This function accounts
%   for atmospheric moisture by using virtual temperature, which corrects
%   for the fact that water vapor is lighter than dry air and thus affects
%   the density and thickness of atmospheric layers.
%
% SYNTAX:
%   Z = computeGeopotentialHeight(T, q, p)
%   Z = computeGeopotentialHeight(T, q, p, Z_sfc)
%
% INPUTS:
%   T     - Temperature [K]
%           Dimensions: (..., nlevels) where last dimension is vertical
%   q     - Specific humidity [kg/kg]
%           Dimensions: must match T
%   p     - Pressure at each level [Pa]
%           Can be either:
%             - 1D array of length nlevels (same pressures everywhere)
%             - Same dimensions as T (pressure varies spatially)
%   Z_sfc - (optional) Surface geopotential height [m]
%           Default is 0 m (sea level)
%           Can be scalar or array matching horizontal dimensions of T
%
% OUTPUTS:
%   Z - Geopotential height [m]
%       Same dimensions as input T
%
% SCIENTIFIC BACKGROUND:
%
%   THE HYPSOMETRIC EQUATION:
%   The hypsometric equation relates the thickness of an atmospheric layer
%   to the mean temperature of that layer. It is derived from the
%   hydrostatic equation and the ideal gas law:
%
%   Hydrostatic balance:    dp = -ρ * g * dz
%   Ideal gas law:          p = ρ * R_d * T_v
%
%   Combining these and integrating from pressure level p1 to p2:
%
%       Z₂ - Z₁ = (R_d / g) * T_v_mean * ln(p₁ / p₂)
%
%   where:
%       Z₁, Z₂ = geopotential heights at pressure levels p₁ and p₂
%       R_d = 287.05 J/(kg·K) = specific gas constant for dry air
%       g = 9.80665 m/s² = gravitational acceleration
%       T_v_mean = mean virtual temperature of the layer
%       ln(p₁/p₂) = natural logarithm of pressure ratio
%
%   This equation tells us that thicker atmospheric layers (larger Z₂ - Z₁)
%   occur when the layer is warmer (larger T_v_mean), which makes physical
%   sense because warm air is less dense and expands vertically.
%
%   ACCOUNTING FOR MOISTURE - VIRTUAL TEMPERATURE:
%   Water vapor molecules (H₂O, molecular weight ≈ 18 g/mol) are lighter
%   than the average dry air molecule (≈ 29 g/mol). This means moist air
%   is less dense than dry air at the same temperature and pressure.
%
%   Virtual temperature (T_v) is the temperature that dry air would need
%   to have the same density as the moist air:
%
%       T_v = T / (1 - (1 - ε) * q)
%       T_v ≈ T * (1 + 0.608 * q)  [for small q]
%
%   where:
%       ε = R_d / R_v ≈ 0.622 (ratio of gas constants)
%       q = specific humidity [kg/kg]
%
%   For typical atmospheric moisture (q ≈ 0.01 kg/kg), this adds about
%   0.6% to the temperature, which translates to several meters difference
%   in geopotential height. For very moist tropical atmospheres
%   (q ≈ 0.02 kg/kg), the correction can be over 1% or ~100 m at 500 hPa.
%
%   WHY THIS MATTERS:
%   Neglecting the moisture correction would cause systematic errors in
%   height calculations, particularly in:
%   - Tropical regions (high moisture content)
%   - Lower troposphere (where absolute humidity is highest)
%   - Summer months (warmer air holds more moisture)
%
%   INTEGRATION APPROACH:
%   This function integrates upward from the surface (or lowest level)
%   using a layer-by-layer approach:
%
%   1. Calculate virtual temperature at each level
%   2. For each layer, compute mean virtual temperature:
%      T_v_mean = (T_v(k) + T_v(k+1)) / 2
%   3. Apply hypsometric equation to get layer thickness
%   4. Add thickness to height below to get height at upper level
%
%   PRESSURE LEVEL ORDERING:
%   The function assumes pressure levels are ordered from high to low
%   pressure (surface to top of atmosphere). If your ERA5 data has
%   levels ordered differently, you'll need to flip the arrays first.
%
% ALGORITHM:
%   Starting from the surface (or lowest level):
%   for each layer from k to k+1:
%       1. Compute virtual temperatures: T_v(k) and T_v(k+1)
%       2. Average for layer: T_v_mean = (T_v(k) + T_v(k+1)) / 2
%       3. Calculate thickness: ΔZ = (R_d/g) * T_v_mean * ln(p(k)/p(k+1))
%       4. Update height: Z(k+1) = Z(k) + ΔZ
%
% EXAMPLES:
%   Example 1: Simple profile (1D vertical)
%       T = [288; 281; 274; 267; 260];  % Temperature [K]
%       q = [0.010; 0.008; 0.005; 0.003; 0.001];  % Specific humidity
%       p = [100000; 92500; 85000; 70000; 50000];  % Pressure [Pa]
%       Z = computeGeopotentialHeight(T, q, p);
%
%   Example 2: ERA5 data (3D: lon × lat × level)
%       T = ncread('era5.nc', 't');      % [K]
%       q = ncread('era5.nc', 'q');      % [kg/kg]
%       p_levels = ncread('era5.nc', 'level') * 100;  % Convert hPa to Pa
%       Z_surface = ncread('era5.nc', 'z') / 9.80665;  % Convert to m
%       Z = computeGeopotentialHeight(T, q, p_levels, Z_surface);
%
%   Example 3: Compare with/without moisture correction
%       Z_moist = computeGeopotentialHeight(T, q, p);
%       q_dry = zeros(size(q));  % No moisture
%       Z_dry = computeGeopotentialHeight(T, q_dry, p);
%       height_difference = Z_moist - Z_dry;
%
% NOTES:
%   - Pressure levels should be ordered from high to low (surface to TOA)
%   - ERA5 typically provides pressure in hPa; multiply by 100 for Pa
%   - ERA5 provides geopotential (Φ) in m²/s²; divide by g to get height
%   - Standard gravity (g = 9.80665 m/s²) is used, not latitude-dependent
%   - For very accurate work, consider latitude-varying gravity
%   - This function uses layer-mean temperatures; for thick layers,
%     this is less accurate than using geometric mean temperatures
%
% REFERENCES:
%   - Wallace, J. M., & Hobbs, P. V. (2006). Atmospheric Science: An
%     Introductory Survey (2nd ed.). Academic Press. Chapter 3.
%   - Holton, J. R. (2004). An Introduction to Dynamic Meteorology
%     (4th ed.). Academic Press. Chapter 1.
%   - Dutton, J. A. (1976). The Ceaseless Wind: An Introduction to the
%     Theory of Atmospheric Motion. McGraw-Hill. Chapter 2.
%
%


% By Andrew John Buggee

%%

function Z = computeGeopotentialHeight(T, q, p, Z_sfc)


% Input validation
if nargin < 3
    error('At least three inputs required: T, q, and p');
end

if nargin < 4
    Z_sfc = 0;  % Default to sea level
end

% Check that T and q have same dimensions
if ~isequal(size(T), size(q))
    error('Temperature and specific humidity must have same dimensions');
end

% Physical constants
R_d = 287.05;        % Specific gas constant for dry air [J/(kg·K)]
g = 9.80665;         % Standard gravity [m/s²]
epsilon = 0.622;     % Ratio R_d/R_v [dimensionless]

% Get dimensions
sz = size(T);
nlevels = sz(1);   % Number of vertical levels (last dimension)

% Handle pressure input
% if isvector(p) && length(p) == nlevels
%     % p is 1D array - same pressure levels everywhere
%     % Need to expand to match T dimensions for vectorized operations
%     p_expanded = reshape(p, [ones(1, ndims(T)-1), nlevels]);
%     p_expanded = repmat(p_expanded, [sz(1:end-1), 1]);
%     p = p_expanded;
% elseif ~isequal(size(p), size(T))
%     error('Pressure must be either 1D array of length nlevels or same size as T');
% end

if ~isequal(size(p), size(T))
    error('Pressure must be either 1D array of length nlevels or same size as T');
end

% Check pressure ordering (should decrease with index)
if ndims(p) == 1
    if p(1) < p(end)
        error('Pressure levels should be ordered from high to low (surface to TOA)');
    end
end

% Calculate virtual temperature at all levels
% T_v = T / (1 - (1 - epsilon) * q)
T_v = T ./ (1 - (1 - epsilon) .* q);

% Initialize geopotential height array
Z = zeros(size(T));


% Set surface/lowest level height
Z(1) = Z_sfc;  % For scalar Z_sfc with vector inputs


% Integrate upward through layers using hypsometric equation
for k = 1:(nlevels-1)
    
    % Extract T_v and p for lower and upper levels (simple vector indexing)
    T_v_lower = T_v(k);
    T_v_upper = T_v(k+1);
    p_lower = p(k);
    p_upper = p(k+1);

    % Calculate mean virtual temperature for the layer
    % Using arithmetic mean (adequate for typical layer thicknesses)
    T_v_mean = (T_v_lower + T_v_upper) / 2;

    % Apply hypsometric equation to get layer thickness
    % ΔZ = (R_d / g) * T_v_mean * ln(p_lower / p_upper)
    thickness = (R_d / g) * T_v_mean * log(p_lower / p_upper);   % meters

    % Add thickness to height at lower level to get height at upper level
    Z(k+1) = Z(k) + thickness;   % meters
end




end