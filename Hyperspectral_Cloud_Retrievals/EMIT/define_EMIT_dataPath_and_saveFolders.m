%%

% (1) folder_extension_number - a number that is used to create a new
% folder using the same nomenclature with 'folder_extension_number'
% appended at the end of the folder. Input of 0 uses the base folder name

% By Andrew John Buggee

%%

function [folder_paths] = define_EMIT_dataPath_and_saveFolders(folder_extension_number)

% Determine which computer you're using
folder_paths.which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.


if folder_extension_number==0



    if strcmp(whatComputer,'anbu8374')==true

        % ------ Folders on my Mac Desktop --------

        % Define the EMIT data folder path

        folder_paths.emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/EMIT/';


        % Define the libRadtran data files path. All paths must be absolute in
        % the INP files for libRadtran
        folder_paths.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';

        % water cloud file location
        folder_paths.water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';

        % mie folder location
        folder_paths.libRadtran_mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';


        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


    elseif strcmp(whatComputer,'andrewbuggee')==true

        % ------ Folders on my Macbook --------

        % Define the EMIT data folder path

        folder_paths.emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/EMIT/INP_OUT/';

        % mie folder location
        folder_paths.libRadtran_mie_folder = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';


        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/';

        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc/';

        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';




        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


    elseif strcmp(whatComputer,'curc')==true


        % ------ Folders on the CU Supercomputer /projects folder --------


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/EMIT/INP_OUT/';


        % Define the EMIT data folder path
        folder_paths.emitDataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


        % libRadtran data folder
        folder_paths.libRadtran_data = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod/';


        % mie folder location
        folder_paths.libRadtran_mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/';


    end





else


    % ----------------------------------------------------------------------
    % --- Append the number at the end and create a new folder if needed ---
    % ----------------------------------------------------------------------
    if strcmp(whatComputer,'anbu8374')==true

        % ------ Folders on my Mac Desktop --------

        % Define the EMIT data folder path

        folder_paths.emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


        % Define the libRadtran data files path. All paths must be absolute in
        % the INP files for libRadtran
        folder_paths.libRadtran_data = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/EMIT/',...
            num2str(folder_extension_number), '/'];

        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/emit_',...
            num2str(folder_extension_number), '/'];


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc_',...
            num2str(folder_extension_number), '/'];


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod_',...
            num2str(folder_extension_number), '/'];


        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


    elseif strcmp(whatComputer,'andrewbuggee')==true

        % ------ Folders on my Macbook --------

        % Define the EMIT data folder path

        folder_paths.emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';

        % libRadtran data folder
        folder_paths.libRadtran_data = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/';


        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/EMIT/INP_OUT/',...
            num2str(folder_extension_number), '/'];


        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/emit_',...
            num2str(folder_extension_number), '/'];


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc_',...
            num2str(folder_extension_number), '/'];

        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod_',...
            num2str(folder_extension_number), '/'];




        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


    elseif strcmp(whatComputer,'curc')==true


        % ------ Folders on the CU Supercomputer /projects folder --------




        % Define the EMIT data folder path
        folder_paths.emitDataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


        % Define the folder path where the mat files of reflectances will be
        % saved
        folder_paths.reflectance_calcs = folder_paths.emitDataPath;


        % libRadtran data folder
        folder_paths.libRadtran_data = '/projects/anbu8374/software/libRadtran-2.0.5/data/';

        % Define the folder path where all .INP files will be saved
        folder_paths.libRadtran_inp = ['/scratch/alpine/anbu8374/EMIT/INP_OUT_',...
            num2str(folder_extension_number), '/'];


        % water cloud file location
        folder_paths.libRadtran_water_cloud_files = ['/projects/anbu8374/software/libRadtran-2.0.5/data/wc_',...
            num2str(folder_extension_number), '/'];


        % libRadtran atmosphere folder location
        folder_paths.atm_folder_path = ['/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod_',...
            num2str(folder_extension_number), '/'];


        % mie folder location
        folder_paths.libRadtran_mie_folder = ['/scratch/alpine/anbu8374/Mie_Calculations/emit_',...
            num2str(folder_extension_number), '/'];


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
