%% run different mie calculations by varrying the wavelength and the alpha parameter

% Each time, generate a netCDF file for use within libRadtran

% By Andrew John Buggee
%%

clear variables

wl = 300:10:2500;    % 221 different wavelengths

% from the VOCALS-REx in-situ measurements, we looked at all alpha
% parameters measured for 73 in-situ measured droplet profiles. The mean
% value was 26.7, the median was 18.9, and the mode was 11.1. We exluced
% nan's in this calculation. 

alpha_param = [1:30, 35:5:60, 80, 100, 125, 150];   % 40 different alpha values

r_eff = 1:1:25;           % microns
% r_eff = 1:1:5;           % microns


%% create mie files

% Each monochromatic mie file takes about 20 minutes to run. I could
% request 221 unique jobs so that all wavelengths run simultaneously. Each
% job would step through the different alpha values. 

% From there, I need to combine netCDF files together so that a single file
% spans the full wavelength and effective radius range for a single alpha
% parameter. In the end, I should have a set of netCDF files equal to the
% length of the alpha parameter.

% Once I have this set of tables, I would need to select the one with the
% alpha value closest to the value I'd like to use. 

% So, let's start by creating 221 directories, one for each monochromatic
% wavelength. Wihtin each directiry will be a set of 40 mie files to run,
% one for each alpha parameter



% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------
    
    folder_path = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Radiative_Transfer_Physics/mieTables_gamma/netCDF_gammaDist/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    folder_path = ['/projects/anbu8374/Matlab-Research/Radiative_Transfer_Physics/',...
        'mieTables_gamma/netCDF_gammaDist/'];

end


% set the mie file parameters
% Set the mie file parameters
mie_program = 'MIEV0';               % type of mie algorithm to run
index_of_refraction = 'water';    % This function only deals with liquid water clouds
size_distribution = 'lognormal';     % Example size distribution type
mie_radius = [r_eff(1), r_eff(end), 1];    % microns
err_msg_str = [];
create_netCDF = true;

for ww = 1:length(wl)

    wl_folder_path = [folder_path, 'wavelength_' num2str(wl(ww)), 'nm'];
    mkdir(wl_folder_path);  % Create directory for each wavelength

    % define the wavelength
    lambda = wl(ww);        % nanometers

    for aa = 1:length(alpha_param)
        % Generate Mie file command for each alpha parameter
        % Create a mie file

        % define the size distribution and the alpha parameter
        size_distribution = {'gamma', alpha_param(aa)};           % droplet distribution

        [input_filename, output_filename] = write_mie_file(mie_program, index_of_refraction,...
            mie_radius, lambda, size_distribution, err_msg_str, aa, alpha_param(aa),...
            [wl_folder_path, '/'], create_netCDF);


    end

end
