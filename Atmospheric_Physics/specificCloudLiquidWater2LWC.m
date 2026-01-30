% specificCloudLiquidWater2LWC - Convert specific cloud liquid water to LWC
%
% DESCRIPTION:
%   Converts ERA5 specific cloud liquid water content to liquid water
%   content (LWC). Specific cloud liquid water is the mass of cloud liquid
%   water droplets per kilogram of total moist air, while LWC is the mass
%   of liquid water per unit volume of air.
%
% SYNTAX:
%   LWC = specificCloudLiquidWater2LWC(q_l, T, q, p)
%
% INPUTS:
%   q_l - Specific cloud liquid water content [kg/kg]
%         ERA5 variable: 'clwc' (cloud liquid water content)
%   T   - Temperature [K]
%         ERA5 variable: 't'
%   q   - Specific humidity [kg/kg]
%         ERA5 variable: 'q'
%         Needed to calculate virtual temperature
%   p   - Pressure [Pa]
%         ERA5 levels in hPa should be multiplied by 100
%
% OUTPUTS:
%   LWC - Liquid Water Content [kg/m³]
%         To convert to g/m³, multiply by 1000
%
% SCIENTIFIC BACKGROUND:
%
%   DEFINITION OF SPECIFIC CLOUD LIQUID WATER CONTENT (q_l):
%   ERA5 defines this as the ratio of liquid water mass to total moist
%   air mass:
%       q_l = m_liquid / m_total
%
%   where m_total includes: dry air + water vapor + cloud liquid +
%   cloud ice + rain + falling snow
%
%   DEFINITION OF LIQUID WATER CONTENT (LWC):
%   LWC is the mass of liquid water per unit volume of air:
%       LWC = m_liquid / V
%
%   CONVERSION RELATIONSHIP:
%   The relationship between these quantities involves the moist air
%   density (ρ_air):
%       LWC = q_l × ρ_air
%
%   This is analogous to converting specific humidity to water vapor
%   density, but for liquid water instead of vapor.
%
%   CALCULATING MOIST AIR DENSITY:
%   The moist air density is calculated using the ideal gas law with
%   virtual temperature:
%       ρ_air = p / (R_d × T_v)
%
%   where:
%       R_d = 287.05 J/(kg·K) is the gas constant for dry air
%       T_v is the virtual temperature [K]
%
%   VIRTUAL TEMPERATURE CONSIDERATIONS:
%   Strictly speaking, virtual temperature should account for all
%   condensed phases (liquid, ice, precipitation). However:
%   - Water vapor typically: q ~ 0.001-0.020 kg/kg
%   - Cloud liquid water typically: q_l ~ 0.0001-0.001 kg/kg
%   - Cloud ice typically: q_i ~ 0.00001-0.0005 kg/kg
%
%   Since condensed water mixing ratios are typically 10-100 times
%   smaller than water vapor mixing ratio, using the standard virtual
%   temperature (based on water vapor only) introduces negligible error
%   (<0.1%) in most cases:
%       T_v = T / (1 - (1 - ε) × q)
%
%   where ε = 0.622 (ratio of gas constants)
%
%   TYPICAL VALUES:
%   - Non-precipitating stratocumulus: LWC ~ 0.1-0.5 g/m³
%   - Cumulus clouds: LWC ~ 0.2-1.0 g/m³
%   - Thick stratus: LWC ~ 0.3-0.6 g/m³
%   - Cumulonimbus: LWC ~ 1-3 g/m³ (can be higher in cores)
%
%   PHYSICAL INTERPRETATION:
%   LWC directly tells you how much liquid water is suspended in a given
%   volume of air. This is crucial for:
%   - Cloud optical properties (extinction, scattering)
%   - Precipitation formation processes
%   - Aircraft icing hazards
%   - Radiative transfer calculations
%
% ALGORITHM:
%   1. Calculate virtual temperature from T and q
%   2. Calculate moist air density: ρ_air = p / (R_d × T_v)
%   3. Calculate LWC: LWC = q_l × ρ_air
%
% EXAMPLES:
%   Example 1: Single profile calculation
%       q_l = 0.0003;    % 0.3 g/kg specific cloud liquid water
%       T = 280;         % 280 K temperature
%       q = 0.005;       % 5 g/kg specific humidity
%       p = 85000;       % 850 hPa pressure
%       LWC = specificCloudLiquidWater2LWC(q_l, T, q, p);
%       % Result: LWC ≈ 0.000324 kg/m³ or 0.324 g/m³
%
%   Example 2: ERA5 data arrays
%       q_l = ncread('era5.nc', 'clwc');  % [kg/kg]
%       T = ncread('era5.nc', 't');       % [K]
%       q = ncread('era5.nc', 'q');       % [kg/kg]
%       p_levels = ncread('era5.nc', 'level') * 100;  % Convert hPa to Pa
%       LWC = specificCloudLiquidWater2LWC(q_l, T, q, p_levels);
%       LWC_gm3 = LWC * 1000;  % Convert to g/m³
%
%   Example 3: Find cloudy regions
%       LWC = specificCloudLiquidWater2LWC(q_l, T, q, p);
%       cloudy_mask = LWC > 0.0001;  % Threshold of 0.1 g/m³
%       mean_LWC_in_clouds = mean(LWC(cloudy_mask));
%
% NOTES:
%   - ERA5 q_l can be very small or zero in cloud-free regions
%   - Typical cloud liquid water: 0.1-0.001 kg/kg (0.1-1.0 g/kg)
%   - For ice clouds, use 'ciwc' (cloud ice water content) instead
%   - Total condensed water = cloud liquid + cloud ice
%   - Pressure must be in Pa; multiply ERA5 levels (hPa) by 100
%   - Temperature must be in Kelvin
%
% RELATED CONVERSIONS:
%   For other hydrometeor types in ERA5:
%   - Cloud ice: Use 'ciwc' with same conversion approach
%   - Rain: Use 'crwc' (rain water content)
%   - Snow: Use 'cswc' (snow water content)
%
% REFERENCES:
%   - Wallace, J. M., & Hobbs, P. V. (2006). Atmospheric Science: An
%     Introductory Survey (2nd ed.). Academic Press. Chapter 6.
%   - Rogers, R. R., & Yau, M. K. (1989). A Short Course in Cloud Physics
%     (3rd ed.). Pergamon Press. Chapter 2.
%   - ERA5 Documentation: https://confluence.ecmwf.int/display/CKB/ERA5
%
% See also: specificHumidity2waterVaporDensity, 
%           specificCloudIceWater2IWC
%
% By Andrew John Buggee




%%

function LWC = specificCloudLiquidWater2LWC(q_l, T, q, p)


% Input validation
if nargin < 4
    error('Four inputs required: q_l, T, q, and p');
end

% Check dimensions match
if ~isequal(size(q_l), size(T)) || ~isequal(size(q_l), size(q))
    error('q_l, T, and q must have the same dimensions');
end

% Handle pressure input (can be scalar, vector, or array matching others)
if isscalar(p)
    p = repmat(p, size(q_l));
elseif isvector(p) && length(p) == size(q_l, ndims(q_l))
    % Expand pressure vector to match dimensions
    sz = size(q_l);
    p_expanded = reshape(p, [ones(1, ndims(q_l)-1), length(p)]);
    p = repmat(p_expanded, [sz(1:end-1), 1]);
elseif ~isequal(size(p), size(q_l))
    error('Pressure dimensions must match q_l or be a compatible vector');
end

% Physical constants
R_d = 287.05;        % Specific gas constant for dry air [J/(kg·K)]
epsilon = 0.622;     % Ratio R_d/R_v [dimensionless]

% Calculate virtual temperature
% T_v = T / (1 - (1 - epsilon) * q)
T_v = T ./ (1 - (1 - epsilon) .* q);

% Calculate moist air density using ideal gas law
% ρ_air = p / (R_d × T_v)
rho_air = p ./ (R_d .* T_v);  % [kg/m³]

% Calculate liquid water content
% LWC = q_l × ρ_air
LWC = q_l .* rho_air;  % [kg/m³]

% Ensure output has same dimensions as input
LWC = reshape(LWC, size(q_l));

end