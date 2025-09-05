%% Generate many measurements with different optical depths, cloud droplet size at top and bottom, and different total column water vapor amounts


clear variables


%% Define the path location where INP files will be stored, and where Reflectances will be stored

clear inputs

% Determine which computer this is being run on
inputs.which_computer = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(inputs.which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];


    % Define the folder path where all .INP files will be saved
    inputs.libRadtran_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';


    % water cloud file location
    inputs.water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';



elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];


    % Define the folder path where all .INP files will be saved
    inputs.libRadtran_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];


    % water cloud file location
    inputs.water_cloud_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/'];

    % mie folder location
    inputs.libRadtran_mie_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/',...
        'libRadtran-2.0.4/Mie_Calculations/'];





elseif strcmp(inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
        'Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/'];



    % Define the folder path where all .INP files will be saved
    inputs.libRadtran_inp = '/scratch/alpine/anbu8374/hyperspectral_retrieval/HySICS/INP-OUT/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';

    % water cloud file location
    inputs.water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';

    % mie folder location
    inputs.libRadtran_mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/';


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

            % Leave a core for overhead
            parpool(p.NumWorkers -1);

        elseif p.NumWorkers<=10

            % Leave a core for overhead
            parpool(p.NumWorkers -1);


        end

    end



end


% If the INP folder path doesn't exist, create a new directory
if ~exist(inputs.libRadtran_inp, 'dir')

    mkdir(inputs.libRadtran_inp)

end


% If the reflectances folder path doesn't exist, create a new directory
if ~exist(inputs.folderpath_reflectance, 'dir')

    mkdir(inputs.folderpath_reflectance)

end


% If the water cloud folder path doesn't exist, create a new directory
if ~exist(inputs.water_cloud_folder_path, 'dir')

    mkdir(inputs.water_cloud_folder_path)

end


%%  Delete old files?
% First, delete files in the HySICS folder
% delete([inputs.libRadtran_inp, '*.INP'])
% delete([inputs.libRadtran_inp, '*.OUT'])

% delete old wc files
% delete([inputs.water_cloud_folder_path, '*.DAT'])

%%

% define the range of independent parameters



r_top = 10;
r_bot = 5;
tau_c = [5,11,17,23];
tcpw = [8, 14, 20];





% ----- unpack parallel for loop variables ------
% We want to avoid large broadcast variables!
libRadtran_inp = inputs.libRadtran_inp;
libRadtran_data_path = inputs.libRadtran_data_path;
wc_folder_path = inputs.water_cloud_folder_path;
which_computer = inputs.which_computer;




% Do you want to model 2 parameters for the droplet profile (r_top and r_bot)
% or 3 (r_top, r_middle, r_bot)

inputs.RT.num_re_parameters = 2;

% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file
print_libRadtran_err = false;

% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';

% we model two free parameters, r_top and r_bot
inputs.RT.num_re_parameters = 2;

[inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, false, [], 'exact', print_libRadtran_err);

%% NO ERROR FILES!

inputs.RT.errMsg = 'quiet';

%% Define the geometry

% Numbers based off NOAA solar calculator (27 Jan 2024 off coast of Chile)
inputs.RT.sza = 10;         % degrees - Solar Zenith Angle
inputs.RT.phi0 = 91.45;        % degrees - Solar Azimuth Angle
inputs.RT.vza = 7;        % degrees - Viewing Zenith Angle
inputs.RT.vaz = 210;       % degrees - Viewing Azimuth Angle




%% Set the total column water vapor?

inputs.RT.modify_total_columnWaterVapor = true;             % modify the full column

inputs.RT.modify_aboveCloud_columnWaterVapor = false;         % don't modify the column above the cloud



%%

% num wavelengths
num_wl = length(inputs.bands2run);

% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);


num_INP_files = num_rTop * num_rBot * num_tauC * num_tcpw * num_wl;

inputFileName = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);


inputs.calc_type = 'simulated_spectra';




% changing variable steps through rTop, rBot, tauC, tcpw, and wavelength
% in for loop speak, it would be:
% for rt = 1:num_rTop
%   for rb = 1:num_rBot
%       for tt = 1:num_tauC
%           for pw = 1:num_tcpw
%               for ww = 1:num_wl
changing_variables_allStateVectors = [reshape(repmat(r_top, num_rBot * num_tauC * num_tcpw * num_wl,1), [],1),...
    repmat(reshape(repmat(r_bot, num_tauC * num_tcpw * num_wl,1), [],1), num_rTop, 1),...
    repmat(reshape(repmat(tau_c, num_tcpw * num_wl,1), [],1), num_rBot * num_rTop, 1),...
    repmat(reshape(repmat(tcpw,  num_wl,1), [],1), num_rBot * num_rTop * num_tauC, 1),...
    repmat(inputs.RT.wavelengths2run, num_rTop * num_rBot * num_tauC * num_tcpw, 1)];


% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables_allStateVectors = [changing_variables_allStateVectors, repmat((1:num_wl)', num_rTop * num_rBot * num_tauC * num_tcpw, 1)];

% First, write all the wc files
temp_names = cell(num_rTop * num_rBot * num_tauC, 1);
wc_filename = cell(num_INP_files, 1);





% only jump on indexes where there is a unique r_top, r_bot and tau pair
idx_unique_wcFiles_idx = 1:(num_wl * num_tcpw):num_INP_files;

parfor nn = 1:length(idx_unique_wcFiles_idx)
    % for nn = 1:length(idx_unique_wcFiles_idx)
    % --------------------------------------
    % ---- Write all Water Cloud files! ----
    % --------------------------------------

    re = create_droplet_profile2([changing_variables_allStateVectors(idx_unique_wcFiles_idx(nn), 1),...
        changing_variables_allStateVectors(idx_unique_wcFiles_idx(nn), 2)],...
        inputs.RT.z, inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector


    temp = write_wc_file(re, changing_variables_allStateVectors(idx_unique_wcFiles_idx(nn), 3),...
        inputs.RT.z_topBottom,inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
        inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
        inputs.RT.indVar, inputs.compute_weighting_functions, inputs.which_computer,...
        idx_unique_wcFiles_idx(nn), inputs.RT.num_re_parameters, wc_folder_path,...
        inputs.libRadtran_mie_folder);

    temp_names{nn} = temp{1};

end

% the wc_filenames should be the same for different wavelengths and total
% precipitable water
for ww = 0:(num_wl * num_tcpw)-1
    wc_filename(idx_unique_wcFiles_idx + ww) = temp_names;
end





% Now write all the INP files
parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


    % set the wavelengths for each file
    wavelengths = changing_variables_allStateVectors(nn, end-2:end-1);

    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the nn value so
    % they can be traced, and are writen over in memory


    inputFileName{nn} = [num2str(mean(wavelengths)), '_','nm_rTop_', num2str(changing_variables_allStateVectors(nn,1)),...
        '_rBot_', num2str(changing_variables_allStateVectors(nn,2)),'_tauC_', num2str(changing_variables_allStateVectors(nn,3)), '_tcpw_',...
        num2str(changing_variables_allStateVectors(nn,4)),'mm.INP'];



    outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


    % ------------------ Write the INP File --------------------
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName{nn}, inputs,...
        wavelengths, wc_filename{nn}, [], [], [], changing_variables_allStateVectors(nn, 4));


end










%% Calculate Reflectance

% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(inputs.RT.wavelengths2run,[],"all"), max(inputs.RT.wavelengths2run, [], "all")];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = inputs.RT.source_file_resolution/3;   % nm



% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
Refl_model_allStateVectors = zeros(num_INP_files, 1);


parfor nn = 1:num_INP_files


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------

    % compute INP file
    runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn},which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(changing_variables_allStateVectors(nn,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables_allStateVectors(nn, end-1) + wl_perturb);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    [Refl_model_allStateVectors(nn), ~] = reflectanceFunction_ver2(inputs, ds,...
        source_flux(idx_wl), spec_response_value(changing_variables_allStateVectors(nn,end),:));



end


toc


%% Rearrange the reflectances

Refl_model_allStateVectors = reshape(Refl_model_allStateVectors, num_wl, []);



%% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
% inputs.measurement.uncert = [0.003, 0.01:0.01:0.1];
inputs.measurement.uncert = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between
% +/- 5% of the measured value

% To sample a normal distribution with mean mu, and standard deviation s,
% we compute the following: y = s * randn() + mu

% We define the standard deviation as the measurement uncertainty divided
% by three. Therefore, after sample a large number of times, 99% of
% measurements will be within +/- measurement uncertainy of the mean

if any(inputs.measurement.uncert > 0)

    inputs.measurement.standard_dev = inputs.measurement.uncert/3;       % this is still just a fraction

    for uu = 1:length(inputs.measurement.uncert)

        clear Refl_model_with_noise_allStateVectors Refl_model_uncert_allStateVectors


        Refl_model_with_noise_allStateVectors = (inputs.measurement.standard_dev(uu) .* Refl_model_allStateVectors) .* randn(size(Refl_model_allStateVectors))...
            + Refl_model_allStateVectors;


        % define the synthetic relfectance uncertainty
        Refl_model_uncert_allStateVectors = inputs.measurement.uncert(uu) .* Refl_model_with_noise_allStateVectors;    % 1/sr



    end

end



%%
% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

% Save the version without an measurement uncertainty. Then we can add
% uncertainty and save the new file

if strcmp(inputs.which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    inputs.folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/'];



elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/paper2_variableSweep/'];


elseif strcmp(inputs.which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    inputs.folderpath_2save = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/',...
        Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset'];



end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_2save, 'dir')

    mkdir(inputs.folderpath_2save)

end


% *** save each wavelength grouping as a standalone measurement ***

for nn = 1:(num_INP_files/num_wl)

    % Grab the reflectance group for the given state vector
    Refl_model = Refl_model_allStateVectors(:,nn);
    Refl_model_with_noise = Refl_model_with_noise_allStateVectors(:,nn);
    Refl_model_uncert = Refl_model_uncert_allStateVectors(:,nn);

    % grab the state vector
    changing_variables = changing_variables_allStateVectors((nn*num_wl - (num_wl-1)) : (nn*num_wl) ,:);


    % set the inputs to have the proper state variables
    inputs.RT.r_top = changing_variables(1,1);
    inputs.RT.r_bot = changing_variables(1,2);
    inputs.RT.tau_c = changing_variables(1,3);
    inputs.RT.waterVapor_column = changing_variables(1,4);


    filename = [inputs.folderpath_2save,'simulated_spectra_HySICS_reflectance_',...
        num2str(numel(inputs.bands2run)), 'bands_',num2str(100*inputs.measurement.uncert), '%_uncert',...
        '_rTop_', num2str(changing_variables(1,1)),...
        '_rBot_', num2str(changing_variables(1,2)), '_tauC_', num2str(changing_variables(1,3)),...
        '_tcwv_', num2str(changing_variables(1,4)),'_vza_', num2str(round(inputs.RT.vza)),...
        '_vaz_', num2str(round(inputs.RT.vaz)), '_sza_', num2str(round(inputs.RT.sza)),...
        '_saz_', num2str(round(inputs.RT.phi0)),...
        '_sim-ran-on-',char(datetime("today")),'.mat'];


    save(filename, "Refl_model", "Refl_model_with_noise", "Refl_model_uncert",...
        "inputs", "spec_response", "changing_variables");

end



