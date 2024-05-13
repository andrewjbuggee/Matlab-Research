%% Compute absorption cross sections from hitran data

% using .par file output from the hitran database, compute the absorption
% cross section [cm^2/molecule]

% ----- INPUTS -----

% (1) T - Temperature of the gas (K)

% (2) P - Pressure of the gas - (atm)

% (3) P_self - partial pressure of the gas in question - (atm)

% (4) wl = wavelengths at which the absorption cross section is desired -
% (nm)

% ----- OUTPUTS -----


% (1) abs_cross_sec - absorption cross section


% By Andrew John Buggee

%%


function [abs_cross_sec] = hitran_compute_abs_cross_section(hitran_file, T, P, P_self, wl)



% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    hitran_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end


%% Check to see if the hitran file is a .par

if strcmp(hitran_file(end-3:end), 'par')==true

    % Then we will import the .par file using the function 'importhitran'
    % and save it as a .mat file

    lines = importhitran([hitran_folder, hitran_file]);
    save([hitran_folder, hitran_file(1:end-3), 'mat'], "lines");

else

    % if its a .mat file, simply load the saved data into the workspave
    load([hitran_folder, hitran_file])

end

%% Specify the wavelength index using the range of interest

wl_index = lines.transitionWavenumber>=(10^4/(wl(end)/1e3)) & lines.transitionWavenumber<=(10^4/(wl(1)/1e3));

%% Define some constants and read the Total Internal Partition Sums

% define the reference temperature
T_ref = 296;            % kelvin

% define the radiometric constants
c2 = 1.4387769;      % cm*K

% Read in the total internal partition Q at the reference temperature for
% the molecule and isotopologue in question
Q = read_total_internal_partition_sum(lines.moleculeNumber(1), lines.isotopologueNumber(1), T);

%% Compute the line intensity

% the line intensity is defined at a reference temperature of 296K
S0 = lines.lineIntensity(wl_index);

% grab the wavenumbers of the desired wavelength range
wavenumber = lines.transitionWavenumber(wl_index);      % cm^-1

% Grab the lower energy energy state of the transition
E_lower = lines.lowerStateEnergy(wl_index);         % cm^-1

% Compute the line intensity at some new temperature
S = S0 .* (Q.ref * exp(-c2 * E_lower./T) .* (1 - exp(-c2*wavenumber./T))) ./...
    (Q.T * exp(-c2 * E_lower./T_ref) .* (1 - exp(-c2*wavenumber./T_ref)));


%% Create a Voigt Lineshape

voigt = voigt_lineShape_for_hitran(lines, wl, T, P, P_self);







end