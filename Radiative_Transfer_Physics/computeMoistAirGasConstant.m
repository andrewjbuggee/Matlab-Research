% Function to compute the moist air gas constant

% INPUTS

%   (1) Specific_humidity


% By Andrew John Buggee


%%



function R_air_moist = computeMoistAirGasConstant(specificHumidity)


% Constants
R_air_dry = 287.052874;               % J/kg/K

% molar mass of dry air
Mol_mass_air_dry = 0.2897;                % kg/mol

% molar mass of water vapor
Mol_mass_h2o_vap = 0.018015;                % kg/mol

% define the epsilon coefficient
eps = Mol_mass_h2o_vap/Mol_mass_air_dry; 




% Calculate the moist air gas constant
R_air_moist = R_air_dry * (1 + ((1-eps)/eps) * specificHumidity);   % J/kg/K


end