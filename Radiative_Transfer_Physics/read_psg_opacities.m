%% This function will read a psg_opc.txt file

% INPUTS:

%   (1) file_name - name of uvspec atm txt file (string) -


% OUTPUTS:
%   (1) cross_section (cm^2/molecule) - a vector containing the altitude levels
%   corresponding to each layer

%   (2) wavelength (nm) - a vector containing pressure at each level




% NOTES:
%   These fiels were created using the planetary spectrum generator hosted
%   by the Goddard Space Flight Center https://psg.gsfc.nasa.gov/

% By Andrew John Buggee
%%

function [wavelength, cross_section] = read_psg_opacities(file_name)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there is 1 input and its the correct file name


if nargin~=1
    error([newline,'Not enough inputs. Need 1: the filename', newline])
end



% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true
    
    psg_folder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/';
    
elseif strcmp(computer_name,'anbu8374')==true
    
    error('You havent stored the atm profiles on you desktop yet!')
    psg_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/';
    
end


%%

% ------------------------------------------------------------
% ---------------------- Read Source File --------------------
% ------------------------------------------------------------


% ---- Open the File ----

file_id = fopen([psg_folder,file_name], 'r');   % 'r' tells the function to open the file for reading

% ------------------------------------------------------
% -------- Reading .txt file using fscanf --------------
% ------------------------------------------------------

% define how the data is written
% '%' starts a new character
% '%*s' tells the code to skip string characters
% 'f' tells us to look for floating point numbers
% 'e' tells the code to look for exponential numbers
% '5.1' tells us to look for a floating point number with 5 entities, 1
% of which comes after the decimal point

%     format_spec = '  %f %e';
%     shape_output = [2 Inf];
%
%     % lets skip the first 11 lines, since they are commented
%     for ii = 1:11
%         fgets(file_id);
%     end
%
%     % now the file pointer will be at the data
%     A = fscanf(file_id, format_spec, shape_output)'; % sxtract data!

% ------------------------------------------------------
% -------- Reading .txt file using textscan ------------
% ------------------------------------------------------
% Or we could use the textscan() function instead, which allows us to define comments to ignore

format_spec = '%f %f';                                  % two floating point numbers
B = textscan(file_id, format_spec, 'Delimiter',' ',...
    'MultipleDelimsAsOne',1, 'CommentStyle','#');


wavelength = B{1};                % wavelength vector

cross_section = B{2};             % molecular cross section vector







end
