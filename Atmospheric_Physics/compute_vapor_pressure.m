% Function to compute the the vapor pressure, the partial pressure of water
% vapor, from the dew point

% INPUTS

%   (1) dew point temperature


% By Andrew John Buggee


%%



function vapor_pressure = compute_vapor_pressure(dew_point_temp)


% Constants
e0 = 6.112; % mb

b = 17.67;

T1 = 273.15;    % K

T2 = 29.65;     % K


% Calculate the vapor pressure
vapor_pressure = e0 .* exp((b * (dew_point_temp - T1)) ./ (dew_point_temp + T2));


end