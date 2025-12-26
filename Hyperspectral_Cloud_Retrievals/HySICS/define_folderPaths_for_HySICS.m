%% Define all the folder paths necessary to read HySICS data, run libRadtran using HySICS data, and save output


% ----- INPUTS ------

% (1) folder_extension_number - a number that is used to create a new
% folder using the same nomenclature with 'folder_extension_number'
% appended at the end of the folder. Input of 0 uses the base folder name


% By Andrew John Buggee

%%

function [folder_paths] = define_folderPaths_for_HySICS(folder_extension_number)



% Determine which computer you're using
folder_paths.which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.


if folder_extension_number==0


    % ---------------------------------------
    % ---- Just use the base folder name ----
    % ---------------------------------------

    if strcmp(folder_paths.which_computer,'anbu8374')==true

        % -----------------------------------------
        % ------ Folders on my Mac Desktop --------
        % -----------------------------------------

        % ***** Define the HySICS Folder with the simulated measurements *****
        folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
            'HySICS/Simulated_spectra/paper2_variableSweep/'];

        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
            'HySICS/Droplet_profile_retrievals/'];

        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/';


        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';

        % mie folder location
        folder_paths.libRadtran_mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';




    elseif strcmp(folder_paths.which_computer,'andrewbuggee')==true



        % -------------------------------------
        % ------ Folders on my Macbook --------
        % -------------------------------------


        % ***** Define the HySICS Folder with the simulated measurements *****
        % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

        folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/',...
            'log_newCov_subset_allBands_VR_inSitu/'];


        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/',...
            'paper2_variableSweep/test_logSpace_newCov_with_VR_inSitu_meas/'];


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/HySICS/INP_OUT/';

        % mie folder location
        folder_paths.libRadtran_mie_folder = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';


        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc/';

        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';





    elseif strcmp(folder_paths.which_computer,'curc')==true


        % ------------------------------------------------
        % ------ Folders on the CU Super Computer --------
        % ------------------------------------------------


        % Define the HySICS simulated spectrum folder

        folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
            'Simulated_spectra/paper2_variableSweep/'];


        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/';


        % libRadtran data folder
        folder_paths.libRadtran_data = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod/';


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/HySICS/INP_OUT/';


        % mie folder location
        folder_paths.libRadtran_mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/';



    end









else


    % ----------------------------------------------------------------------
    % --- Append the number at the end and create a new folder if needed ---
    % ----------------------------------------------------------------------


    if strcmp(folder_paths.which_computer,'anbu8374')==true

        % -----------------------------------------
        % ------ Folders on my Mac Desktop --------
        % -----------------------------------------

        % ***** Define the HySICS Folder with the simulated measurements *****
        folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
            'HySICS/Simulated_spectra/paper2_variableSweep/'];

        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
            'HySICS/Droplet_profile_retrievals/'];

        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


        
        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/INP_OUT_',...
                                        num2str(folder_extension_number), '/'];

        
        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations_',...
                                                num2str(folder_extension_number), '/'];


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod_',...
                                        num2str(folder_extension_number), '/'];


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc_',...
                                                        num2str(folder_extension_number), '/'];

        




    elseif strcmp(folder_paths.which_computer,'andrewbuggee')==true



        % -------------------------------------
        % ------ Folders on my Macbook --------
        % -------------------------------------


        % ***** Define the HySICS Folder with the simulated measurements *****
        % folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        %     'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

        folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/',...
            'log_newCov_subset_allBands_VR_inSitu_1/'];


        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
            'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/',...
            'paper2_variableSweep/test_logSpace_newCov_with_VR_inSitu_meas/'];


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/HySICS/INP_OUT_',...
                                                        num2str(folder_extension_number), '/'];


        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/',...
            'Mie_Calculations_', num2str(folder_extension_number), '/'];



        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc_',...
                                                        num2str(folder_extension_number), '/'];


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod_',...
                                                        num2str(folder_extension_number), '/'];






    elseif strcmp(folder_paths.which_computer,'curc')==true


        % ------------------------------------------------
        % ------ Folders on the CU Super Computer --------
        % ------------------------------------------------


        % Define the HySICS simulated spectrum folder

        folder_paths.HySICS_simulated_spectra = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
            'Simulated_spectra/paper2_variableSweep/'];


        % ---- Define where the retrievals will be stored ---
        folder_paths.HySICS_retrievals = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/';


        % libRadtran data folder
        folder_paths.libRadtran_data = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/projects/anbu8374/software/libRadtran-2.0.5/data/wc_',...
                                                        num2str(folder_extension_number), '/'];



        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod_',...
                                                        num2str(folder_extension_number), '/'];



        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/scratch/alpine/anbu8374/HySICS/INP_OUT_', num2str(folder_extension_number), '/'];



        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/scratch/alpine/anbu8374/Mie_Calculations/Mie_Calculations_', num2str(folder_extension_number), '/'];




    end









end






% If the libRadtran INP folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end


% If the mie INP folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_mie_folder, 'dir')

    mkdir(folder_paths.libRadtran_mie_folder)

end



% If the water cloud file path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_water_cloud_files, 'dir')

    mkdir(folder_paths.libRadtran_water_cloud_files)

end


% If the atmosphere file path doesn't exit, create a new directory
if ~exist(folder_paths.atm_folder_path, 'dir')

    mkdir(folder_paths.atm_folder_path)

end





end

