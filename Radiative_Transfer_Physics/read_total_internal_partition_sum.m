%% Read HITRAN's molparam file to get the Total Internal Partition at the reference temperature

% This is for a single moleucle and a single isotopologue at a single
% temperature

% -------------------------------
% ---------- INPUTS -------------
% -------------------------------
% (1) molecule_number - molecule number according to the hitran database.
% This is store in the hitran line parameter files 

% (2) isotopologue_number - isotopologue number according to the hitran database.
% This is store in the hitran line parameter files 

% (3) T - temperature of the gas (K)


% -------------------------------
% ---------- OUTPUTS -------------
% -------------------------------
% (1) Q_tips - a structure containing the total internal partition sums at
% the hitran reference temperature of 296 K and the user input temperature,
% T



% By Andrew John Buggee

function Q_tips = read_total_internal_partition_sum(molecule_number, isotopologue_number, T)

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    hitran_folder = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end

%% Read in the file

% ---- Open the File ----

file_id = fopen([hitran_folder,num2str(molecule_number), '_', num2str(isotopologue_number),...
    '_Q.txt'], 'r');   % 'r' tells the function to open the file for reading



% ------------------------------------------------------
% -------- Reading .dat file using textscan ------------
% ------------------------------------------------------
% Or we could use the textscan() function instead, which allows us to define comments to ignore

format_spec = '%f %f';                                  % two floating point numbers
ds = textscan(file_id, format_spec);


%% Compute the reference TIPS value and the TIPS value at some other temperature

T_ref = 296;        % K - reference temperature according to Hitran

% Store the reference TIPS value
Q_tips.ref = ds{2}(ds{1}==T_ref);

% Compute the TIPS value at some other temperature, T
Q_tips.T = interp1(ds{1}, ds{2}, T);


end
