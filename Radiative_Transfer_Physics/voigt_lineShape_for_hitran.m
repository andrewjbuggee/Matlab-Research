%% Compute the Voigt spectral line shape using hitran data



% By Andrew John Buggee

%%

function voigt = voigt_lineShape_for_hitran(hitran_lines, wavelength_boundaries, T, P, P_self)


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    hitran_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end


%% Compute the Doppler broadened line shape (gaussian line shape)
% Doppler broadening tends to dominate in the upper atmosphere where
% pressure is low

% we have to compute a line shape for each spectral transition

% We need the molar mass of the isotopologue

molar_mass = read_isotopologue_molar_mass_hitran(hitran_lines);     % g/mol
% convert molar_mass to kg/mol
molar_mass = molar_mass/1000;                           % kg/mol

% load physical constants
con = physical_constants();

% compute the doppler profile half-width at half-max
% use SI units for everything, except leave the spectral dimension in
% wavenumbers with units of inverse centimeters
doppler_hwhm = hitran_lines.transitionWavenumber./con.c .* sqrt((2*con.N_A * con.k_B * T * log(2))/molar_mass);


% compute the Doppler line shape for each line transition
for vv = 1:length(hitran_lines.transitionWavenumber)
    
    % Lines are narrow! Some less than a wavenumber
    wavenumber_vector{vv} = linspace(hitran_lines.transitionWavenumber(vv) - 2.5*doppler_hwhm(vv),...
        hitran_lines.transitionWavenumber(vv) + 2.5*doppler_hwhm(vv), 100);

    f_doppler{vv} = sqrt(log(2)./(pi*doppler_hwhm(vv).^2)) .* ...
        exp(- ((wavenumber_vector{vv} - hitran_lines.transitionWavenumber(vv)).^2 * log(2))./doppler_hwhm(vv)^2);

end



%% Compute the Pressure broadened line shape (Lorentz line shape)
% In the lower atmosphere, pressure broadening dominates

% we have to compute a line shape for each spectral transition

% define the hitran reference temperature
T_ref = 296;        % K 

% compute the doppler profile half-width at half-max
% use SI units for everything, except leave the spectral dimension in
% wavenumbers with units of inverse centimeters
pressure_hwhm = (T_ref/T).^hitran_lines.temperatureDependence .* hitran_lines.airBroadenedWidth .* (P - P_self) +...
    hitran_lines.selfBroadenedWidth .* P_self;


% compute the pressure broadened line shape for each line transition
for vv = 1:length(hitran_lines.transitionWavenumber)

     % Lines are narrow! Some less than a wavenumber
    wavenumber_vector{vv} = linspace(hitran_lines.transitionWavenumber(vv) - 2.5*pressure_hwhm(vv),...
        hitran_lines.transitionWavenumber(vv) + 2.5*pressure_hwhm(vv), 100);

    f_pressure{vv} = pressure_hwhm(vv)./pi .* 1./...
        (pressure_hwhm(vv).^2 + (wavenumber_vector{vv} - (hitran_lines.transitionWavenumber(vv) +...
        hitran_lines.pressureShift(vv) .* P).^2));

end


%% The Voigt Lineshape is the convolution of the doppler broadened and pressure broadened line shapes

for vv = 1:length(hitran_lines.transitionWavenumber)

    voigt.wavenumbers{vv} = wavenumber_vector{vv};
    voigt.shape{vv} = conv(f_doppler{vv}, f_pressure{vv});

end


end