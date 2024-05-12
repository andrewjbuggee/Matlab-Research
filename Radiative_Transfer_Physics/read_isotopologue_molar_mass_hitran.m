%% REad the molar mass of the isotoplogue as define by the hitran data base


% Andrew JOhn Buggee

%%

function molar_mass = read_isotopologue_molar_mass_hitran(hitran_lines)


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    error('You havent stored the atm profiles on you desktop yet!')

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end

%% Read in the file

% ---- Open the File ----

file_id = fopen([hitran_folder, 'molparam_Q_reference.txt'], 'r');   % 'r' tells the function to open the file for reading



% ------------------------------------------------------
% -------- Reading .dat file using textscan ------------
% ------------------------------------------------------

% to read in the data in this .txt file, use a for loop. This will allow
% you to read in the data for each molecule individually

N_molecules = 55;

% create a cell array
ds = cell(N_molecules, 1);

% define the initial number of header lines
n_header_lines = 2;
for nn = 1:N_molecules

    if nn>1
        % define the start of data for the next molecule
        n_header_lines = n_header_lines + size(ds{nn-1}.data, 1) +1;
    end

    ds{nn} = importdata([hitran_folder, 'molparam_Q_reference.txt'], ' ', n_header_lines);

end


%% Grab the molar mass associated with the molecular number and istopolague from the lines data

molar_mass = ds{hitran_lines.moleculeNumber(1)}.data(hitran_lines.isotopologueNumber(1), end);      % g/mol


end
