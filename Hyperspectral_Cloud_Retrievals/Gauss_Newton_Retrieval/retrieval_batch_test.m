%% Retreive vertical profiles in a loop using simulated HySICS reflectance measurements
% The loop can be used to change pixels or to change retrieval settings

% This script uses my own TBLUT algorithm as the apriori values



% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%% Define the HySICS folders for the machine you're using



% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the HySICS Folder with the simulated measurements *****
    folder_paths.HySICS_simulated_spectra = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];

    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/'];

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/';


    % water cloud file location
    folder_paths.water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';



elseif strcmp(which_computer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ***** Define the HySICS Folder with the simulated measurements *****
    folder_paths.HySICS_simulated_spectra = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];

    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/HySICS/'];

    % water cloud file location
    folder_paths.water_cloud_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % Define the HySICS simulated spectrum folder

    folder_paths.HySICS_simulated_spectra = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';


    % ---- Define where the retrievals will be stored ---
    folder_paths.HySICS_retrievals = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/';


    % water cloud file location
    folder_paths.water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';

    % Define the folder path where all .INP files will be saved
    folder_paths.libRadtran_inp = '/scratch/alpine/anbu8374/HySICS/INP_OUT/';


    % *** Start parallel pool ***
    % Is parpool running?
    p = gcp('nocreate');
    if isempty(p)==true

        % first read the local number of workers avilabile.
        p = parcluster('local');
        % start the cluster with the number of workers available
        if p.NumWorkers>64
            % Likely the amilan128c partition with 2.1 GB per core
            % Leave some cores for overhead
            parpool(p.NumWorkers - 8);

        elseif p.NumWorkers<=64 && p.NumWorkers>10

            parpool(p.NumWorkers - 2);

        elseif p.NumWorkers<=10

            parpool(p.NumWorkers);

        end

    end



end

% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end





%%   Delete old files?
% First, delete files in the HySICS folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.water_cloud_folder_path, '*.DAT'])



%% LOAD SIMULATED HYSICS DATA

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_66Bands_20mm-aboveCloud-WV_sim-ran-on-08-Jul-2025_rev1';  % sza = 0, vza = 0



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev1.mat']); % sza = 10, vza = 0

    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-12-May-2025_rev1.mat']); % sza = 0, vza = 0

    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev2.mat']); % sza = 20, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev3.mat']); % sza = 30, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev4.mat']); % sza = 40, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev5.mat']); % sza = 50, vza = 0
    % 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-15-May-2025_rev6.mat']); % sza = 60, vza = 0

    % r_top = 9.5, r_bot = 4, tau_c = 6
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-17-Jun-2025_rev1.mat';

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47 bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_47Bands_20mm-totalColumnWaterVapor_sim-ran-on-07-Jul-2025_rev1';

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 66
    % Bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_66Bands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1';
  

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, ALL bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_allBands_20mm-totalColumnWaterVapor_sim-ran-on-08-Jul-2025_rev1';


    % test file with just 5 wavelengths
    % filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_5wavelength_test_sim-ran-on-10-Jun-2025_rev1.mat';




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 47
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    %filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-07-Jul-2025_rev1.mat';


    % r_top = 9.5, r_bot = 4, tau_c = 6, total_column_waterVapor = 20, 66
    % bands
    % simulated calcs for MODIS obs on fig 3.a for paper 1
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_66Bands_20mm-aboveCloud-WV_sim-ran-on-08-Jul-2025_rev1.mat';



end


simulated_measurements = load([folder_paths.HySICS_simulated_spectra,filename]);


%% Create the name of the file to save all output to

rev = 1;

folder_paths.saveOutput_filename = [folder_paths.HySICS_retrievals,'dropletRetrieval_HySICS_', num2str(numel(simulated_measurements.inputs.bands2run)),...
    'bands_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];



while isfile(filename)
    rev = rev+1;
    if rev<10
        filename = [filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        filename = [filename(1:end-6), num2str(rev),'.mat'];
    end
end


