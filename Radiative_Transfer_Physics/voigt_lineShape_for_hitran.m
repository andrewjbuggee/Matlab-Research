%% Compute the Voigt spectral line shape using hitran data



% By Andrew John Buggee

%%

function voigt = voigt_lineShape_for_hitran(hitran_lines, wavelength_boundaries)


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    error('You havent stored the atm profiles on you desktop yet!')

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end


%% Compute the half-width at half-max of the Doppler broadened line (gaussian line shape)

% We need the molar mass of the isotopologue

molar_mass = read_isotopologue_molar_mass_hitran(hitran_lines);





end