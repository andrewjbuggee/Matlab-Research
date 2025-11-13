%% Retreive vertical profiles using simulated HySICS reflectance measurements
% This is a function that retrieves a vertical droplet profile and the
% integrated column water vapor amount above cloud

% This script retrieves 4 variables: ln(r_top), ln(r_bot), ln(tau_c), and acpw


% INPUTS

% (1) 


% By Andrew John Buggee

%%

function [rho, waterVapor_concentration_cm3] = compute_waterVapor_concentration_from_radiosonde(dew_point_temp, temp_env)



% Convert temperatures from Celsius to Kelvin
temp_env = temp_env + 273.15;                % K
dew_point_temp = dew_point_temp + 273.15;    % K

% % for moist air, a mixture of dry air and water vapor molecules, the idea
% % gas law is:
% % rho = P / (R_dry * T_v)
% % where R_dry is the gas constant for dry air and T_v is the virtual
% % temperature.
% 
% % Constants
% R_air_dry = 287.052874;               % J/kg/K
% 
% % molar mass of dry air
% Mol_mass_air_dry = 0.2897;                % kg/mol
% 
% % molar mass of water vapor
% Mol_mass_h2o_vap = 0.018015;                % kg/mol
% 
% eps = Mol_mass_h2o_vap/Mol_mass_air_dry;


% Compute the water vapor partial pressure from the dew point temperature
vapor_pressure = compute_vapor_pressure(dew_point_temp);              % mb - partial pressure of water vapor per unit volume of air
vapor_pressure_pa = vapor_pressure .* 100;                            % Pa

% using the partial pressure of water vapor, calculate the mass of water
% vapor per unit volume of air

% water vapor gas constant
R_waterVapor = 461.52;               % J/(kg K) 
% Avogadro's constant
N_A = 6.02214076e+23;                % molecules/mol
% universal gas constant
R_uni = 8.31446262;                % J/K/mol


rho = vapor_pressure_pa ./(R_waterVapor * temp_env);              % kg / m^3

waterVapor_concentration = (vapor_pressure_pa * N_A) ./ (R_uni * temp_env);         % num particles / m^3

waterVapor_concentration_cm3 = waterVapor_concentration ./ 1e6;         % num particles / cm^3






end