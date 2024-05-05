%% Read HITRAN's molparam file to get the Total Internal Partition at the reference temperature

% This is for a single moleucle and a single isotopologue

% By Andrew John Buggee

function Q_ref = read_reference_total_internal_partition(molecule_number, isotopologue_number)

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
    
    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading
    

    
    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore
    
    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');



end
