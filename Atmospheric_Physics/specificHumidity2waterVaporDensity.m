%% Function to compute water vapor density from specific humidty

% DESCRIPTION:
%   Converts specific humidity (mass mixing ratio of water vapor) to water
%   vapor density (mass of water vapor per unit volume of air). This
%   conversion requires knowledge of the moist air density, which depends
%   on temperature and pressure through the ideal gas law.
%
% SYNTAX:
%   rho_v = specificHumidity2waterVaporDensity(q, T, p)
%   rho_v = specificHumidity2waterVaporDensity(q, T, p, use_virtual_temp)
%
% INPUTS:
%   q                - Specific humidity [kg/kg or g/kg]
%                      If g/kg, function will detect and convert
%   T                - Temperature [K]
%   p                - Pressure [Pa]
%   use_virtual_temp - (optional) Logical flag to use virtual temperature
%                      correction. Default is false (uses approximation)
%                      true:  Uses virtual temperature (more accurate)
%                      false: Uses standard approximation (~1% error)
%
% OUTPUTS:
%   rho_v - Water vapor density [kg/m³]
%           To convert to g/m³, multiply output by 1000
%
% SCIENTIFIC BACKGROUND:
%   Specific humidity (q) is defined as the ratio of water vapor mass (m_v)
%   to the total mass of moist air (m_moist):
%       q = m_v / m_moist
%
%   Water vapor density (rho_v) is the mass of water vapor per unit volume:
%       rho_v = m_v / V
%
%   The relationship between these quantities involves the moist air density:
%       rho_v = q * rho_air
%
%   where rho_air is the density of moist air, calculated from the ideal
%   gas law for moist air.
%
% IDEAL GAS LAW FOR MOIST AIR:
%   The ideal gas law relates pressure, density, and temperature:
%       p = rho_air * R_d * T_v
%
%   Rearranging for moist air density:
%       rho_air = p / (R_d * T_v)
%
%   where:
%       R_d = 287.05 J/(kg·K) is the specific gas constant for dry air
%       T_v is the virtual temperature [K]
%
% VIRTUAL TEMPERATURE:
%   Virtual temperature is the temperature that dry air would need to have
%   the same density as moist air at the same pressure. Water vapor is
%   lighter than dry air (molecular weight 18 vs 29 g/mol), so moist air
%   is less dense than dry air at the same temperature and pressure.
%
%   The virtual temperature is:
%       T_v = T / (1 - (1 - epsilon) * q)
%
%   where epsilon = R_d / R_v ≈ 0.622 is the ratio of gas constants for
%   dry air and water vapor.
%
%   For typical atmospheric moisture levels (q < 0.02 kg/kg), the virtual
%   temperature correction is small (~1% or less), so the approximation
%   T_v ≈ T is often used for simplicity.
%
% APPROXIMATION VS. EXACT CALCULATION:
%   Approximate (default):  rho_v ≈ (q * p) / (R_d * T)
%   Exact:                  rho_v = q * p / (R_d * T_v)
%
%   The approximation introduces ~1% error for typical atmospheric
%   conditions and is computationally simpler.
%
% EXAMPLES:
%   Example 1: Basic usage with ERA5 data
%       q = 0.008;      % 8 g/kg specific humidity
%       T = 288.15;     % 15°C in Kelvin
%       p = 101325;     % Sea level pressure in Pa
%       rho_v = specificHumidity2waterVaporDensity(q, T, p);
%       % Result: rho_v ≈ 0.00693 kg/m³ or 6.93 g/m³
%
%   Example 2: Using virtual temperature correction
%       rho_v_exact = specificHumidity2waterVaporDensity(q, T, p, true);
%
%   Example 3: Processing ERA5 arrays
%       % Load ERA5 data
%       q_era5 = ncread('era5_file.nc', 'q');    % [kg/kg]
%       T_era5 = ncread('era5_file.nc', 't');    % [K]
%       p_era5 = ncread('era5_file.nc', 'level') * 100;  % Convert hPa to Pa
%       rho_v = specificHumidity2waterVaporDensity(q_era5, T_era5, p_era5);
%
% NOTES:
%   - All inputs must have compatible dimensions (scalar or matching arrays)
%   - Function automatically detects if q is in g/kg (values > 0.1) and
%     converts to kg/kg
%   - Pressure must be in Pascals. If you have hPa, multiply by 100
%   - Temperature must be in Kelvin. To convert from Celsius: T_K = T_C + 273.15
%
% REFERENCES:
%   - Wallace, J. M., & Hobbs, P. V. (2006). Atmospheric Science: An
%     Introductory Survey (2nd ed.). Academic Press.
%   - Rogers, R. R., & Yau, M. K. (1989). A Short Course in Cloud Physics
%     (3rd ed.). Pergamon Press.




% By Andrew John Buggee


%%


function [rho_v, waterVapor_concentration_cm3] = specificHumidity2waterVaporDensity(q, T, p, use_virtual_temp)




% Input validation
if nargin < 3
    error('At least three inputs required: q, T, and p');
end

if nargin < 4
    use_virtual_temp = false;  % Default to approximation
end

% Check that inputs have compatible dimensions
if ~isequal(size(q), size(T)) || ~isequal(size(q), size(p))
    error('Input arrays must have the same dimensions');
end

% Physical constants
R_d = 287.05;        % Specific gas constant for dry air [J/(kg·K)]
R_v = 461.5;         % Specific gas constant for water vapor [J/(kg·K)]
epsilon = R_d / R_v; % Ratio of gas constants ≈ 0.622 [dimensionless]

% Check if specific humidity is in g/kg and convert to kg/kg if needed
if any(q(:) > 0.1)
    warning('Specific humidity values > 0.1 detected. Assuming units are g/kg and converting to kg/kg');
    q = q / 1000;  % Convert g/kg to kg/kg
end

% Calculate water vapor density
if use_virtual_temp
    % Exact method: Use virtual temperature correction
    % T_v = T / (1 - (1 - epsilon) * q)
    T_v = T ./ (1 - (1 - epsilon) .* q);
    
    % Calculate moist air density using virtual temperature
    rho_air = p ./ (R_d .* T_v);  % [kg/m³]
    
    % Calculate water vapor density
    rho_v = q .* rho_air;  % [kg/m³]
else
    % Approximate method: Neglect virtual temperature correction
    % This is accurate to within ~1% for typical atmospheric conditions
    rho_v = (q .* p) ./ (R_d .* T);  % [kg/m³]
end

% Ensure output has same dimensions as input
rho_v = reshape(rho_v, size(q));     % [kg/m³]


% lastly, compute the number concentration of water vapor from the density
con = physical_constants;

Na = con.N_A;                  % #/mol - avogadro's constant
M_wv = con.Mol_mass_h2o_vap;   % kg/mol - molar mass of water vapor

waterVapor_concentration = rho_v * (Na/M_wv);         % num particles / m^3

waterVapor_concentration_cm3 = waterVapor_concentration ./ 1e6;         % num particles / cm^3


end

