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


% (1) Brightness Temperature (Tb) - units in Kelvin (K)


% By Andrew John Buggee

%%


function [outputArg1,outputArg2] = hitran_compute_abs_cross_section(hitran_file, T, P, P_self wl)



% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true

    error('You havent stored the atm profiles on you desktop yet!')

elseif strcmp(computer_name,'anbu8374')==true

    hitran_folder = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/HiTran_data/';


end



%% Read in the hitran output file at the desired wavelengths

% ---- Open the File ----
    
file_id = fopen([hitran_folder, hitran_file], 'r');   % 'r' tells the function to open the file for reading

% ------------------------------------------------------
% -------- Reading .txt file using textscan ------------
% ------------------------------------------------------
% Or we could use the textscan() function instead, which allows us to define comments to ignore

format_spec = ['%1d %1d %12.6f %10.3e %10.3e %5.4f %5.3f ',...
    '%10.4f %4.2f %8.6f %15s %15s %15s %15s %1d %2d %1s %7.1f %7.1f'];
shape_output = [2, 19];

% now the file pointer will be at the data
A = fscanf(file_id, format_spec, shape_output); % sxtract data!


Ds = textscan(file_id, format_spec, 'Delimiter',' ',...
    'MultipleDelimsAsOne',1, 'CommentStyle','#');



end