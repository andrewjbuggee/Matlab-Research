%% Optimal estimation of a vertical droplet profile using EMIT data

% Retrieve a droplet profile from EMIT data
% Use spectral reflectances outside of water vapor bands


% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;

%% Define the EMIT data file to use by defining the data folder

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

% emitDataFolder = '17_Jan_2024_coast/';

% 27 january has overlap with MODIS observations
emitDataFolder = '27_Jan_2024/';

% -------------------------------------


%% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders();


%% Load EMIT Data

[emit,L1B_fileName] = retrieveEMIT_data([folder_paths.emitDataPath, emitDataFolder]);


%% Set INP filename

folder_paths.libRadtran_inp = [folder_paths.libRadtran_inp, 'EMIT_', emitDataFolder(1:end-1), '_',...
    L1B_fileName{1}(27:30),'/'];


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

%% Define the pixels to use for the retrieval


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% pixels2use.row = 1112;
% pixels2use.col = 974;

% 17_Jan_2024_coast - small optical depth
% pixels2use.row = [912, 913];
% pixels2use.col = [929, 929];

% 17_Jan_2024_coast - clear sky over ocean
% pixels2use.row = 1021;
% pixels2use.col = 536;

% 17_Jan_2024_coast - my TBLUT algorithm found an optical depth of 6.57 and
% an effective radius of 10.79
% pixels2use.row = 932;
% pixels2use.col = 960;

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932, 932];
% pixels2use.col = [970, 969];

% 17_Jan_2024_coast - optical depth of 3.2
% pixels2use.row = 932;
% pixels2use.col = 970;

% 17_Jan_2024_coast - optical depth of 3.8
% pixels2use.row = 932;
% pixels2use.col = 969;

% 17_Jan_2024_coast - optical depth of 8.7
% pixels2use.row = 969;
% pixels2use.col = 991;

% 17_Jan_2024_coast - optical depths of 6.6, 11.6
% pixels2use(1).row = 932; pixels2use(1).col = 960;
% pixels2use(2).row = 969; pixels2use(2).col = 984;

% 17_Jan_2024_coast - optical depths of 9.22, 11.6, 14.53 19.8
% pixels2use.row = [969, 969, 969, 969];
% pixels2use.col = [989, 984, 980, 957];

% 17_Jan_2024_coast - optical depths of 8.7, 9.22, 9.68, 10.3, 11.6, 12.54,
% 13.61, 14.53, 16.79, 19.8
% pixels2use.row = [969, 969, 969, 969, 969, 969, 969, 969, 969, 969];
% pixels2use.col = [991, 989, 987, 986, 984, 980, 976, 974, 966, 957];


% 17_Jan_2024_coast - optical depths of 3.2, 6.6, 10.3, 12.54, 14.53, 19.8
% pixels2use.row = [932, 932, 969, 969, 969, 969];
% pixels2use.col = [970, 960, 986, 980, 974, 987];

% 17_Jan_2024_coast - optical depths of 10.3
% pixels2use.row = 969;
% pixels2use.col = 986;


% 17_Jan_2024_coast - optical depths of 0(clear sky), 6.6, 9.22, 11.6, 14.53 19.8
% pixels2use.row = [1021, 932, 969, 969, 969, 969];
% pixels2use.col = [536, 960, 989, 984, 980, 957];



% 27_Jan_2024 - ** Overlap with MODIS **
% ** Time difference bu a couple minutes **
% MODIS retrieved an optical depth of 32.63 and
% an effective radius of 13.27
% modis_pixel_row = 1458;
% modis_pixel_col = 1288;
% pixels2use.row = 1242;
% pixels2use.col = 973;


% 27_Jan_2024 - ** Overlap with MODIS **
% ** Time difference bu a couple minutes **
% MODIS retrieved an optical depth of 12.07 and
% an effective radius of 7.94
% modis_pixel_row = 1481;
% modis_pixel_col = 1285;
pixels2use.row = 1242;
pixels2use.col = 640;


% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, [size(emit.radiance.measurements,1),...
    size(emit.radiance.measurements, 2)]);



%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
GN_inputs = create_gauss_newton_inputs_for_emit_ver2(emitDataFolder, folder_paths, L1B_fileName, emit);
%inputs = create_emit_inputs_hyperspectral_top_middle(emitDataFolder, folder2save, L1B_fileName, emit);

% *** Check Inputs ***

%% Set the wavelenghts!!

% --- Indexes using same 35 as above, in addition to 29 water vapor bands ---
% This set has a total of 64 bands. They are not exactly the same set as
% the 66 HySICS bands used to retrieve column water vapor because the
% HySICS channels are more narrow.
GN_inputs.bands2run = [17, 20, 25, 32, 39, 65, 66, 67, 68, 71, 74, 78, 86, 87, 88, 89, 90,...
    94, 97, 99, 101, 105, 115, 116, 117, 139, 141, 142, 145, 147, 148, 149, 151, 156,...
    157, 158, 172, 175, 176, 187, 189, 190, 210, 212, 213, 214, 215, 216, 217, 218, 219,...
    220, 222, 231, 233, 234, 235, 236, 249, 250, 251, 252, 253, 254]';


%% Override input settings with MODIS derived values

% Values for 27_Jan_2024 - ** pixel [1242, 973] **
% override the cloud top height
GN_inputs.RT.z_topBottom = [0.700, 0.500];    % km

% change inputs that depend on z_topBottom
% Water Cloud depth
GN_inputs.RT.H = GN_inputs.RT.z_topBottom(1) - GN_inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud

GN_inputs.RT.z_edges = linspace(GN_inputs.RT.z_topBottom(2), GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers+1);   % km - the edges of each layer
GN_inputs.RT.z = linspace(GN_inputs.RT.z_topBottom(2), GN_inputs.RT.z_topBottom(1), GN_inputs.RT.n_layers);        % km - altitude above ground vector



%% This retrieval does NOT retrieve column water vapor. What should the forward model assumption be?

GN_inputs.RT.modify_total_columnWaterVapor = false;             % don't modify the full column

% *** Retreive Column Water Vapor! ***
GN_inputs.RT.modify_aboveCloud_columnWaterVapor = true;         % modify the column above the cloud

%% Set output filename

rev = 1;


folder_paths.saveOutput_filename = [folder_paths.emitDataPath, emitDataFolder,'Droplet_profile_retrievals/',...
    num2str(numel(GN_inputs.bands2run)),...
    'bands_ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];




while isfile(folder_paths.saveOutput_filename)
    rev = rev+1;
    if rev<10
        folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        folder_paths.saveOutput_filename = [folder_paths.saveOutput_filename(1:end-6), num2str(rev),'.mat'];
    end
end


%% Define the spectral response function of EMIT for the desired Bands

% create the spectral response functions
[GN_inputs, spec_response] = create_EMIT_specResponse(emit, GN_inputs);


%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
GN_inputs = define_source_for_EMIT(GN_inputs, emit);


%% Convert radiance measurements to TOA reflectance for the desired pixels

emit = convert_EMIT_radiance_2_reflectance(emit, GN_inputs);


%% Compute the radiance measurement uncertainty

[emit.radiance.uncertainty, emit.radiance.uncertainty_percent_perChannel] = compute_EMIT_radiance_uncertainty(emit);


%% Compute the reflectance uncertainty

emit.reflectance.uncertainty = compute_EMIT_reflectance_uncertainty(emit, GN_inputs);



%%  *** Start parallel pool ***

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

        parpool(p.NumWorkers);

    elseif p.NumWorkers<=10

        parpool(p.NumWorkers);

    end

end




%% Check the thermodynamic phase of the defined pixels

%GN_inputs.cloudPhase = determine_cloud_phase_emit(emit, pixels2use);


%% Compute the TBLUT retrieval estimate

tic
tblut_retrieval = TBLUT_forEMIT(emit, spec_response, emitDataFolder, folder_paths);
disp([newline, 'TBLUT retrieval took ', num2str(toc), 'seconds to run', newline])




%% Create the Model and Measurement prior


GN_inputs = create_model_prior_covariance_EMIT_top_bottom_ver2(GN_inputs, tblut_retrieval, true);
%inputs = create_model_prior_covariance_EMIT_top_middle(inputs, pixels2use, tblut_retrieval, true);

GN_inputs = create_EMIT_measurement_covariance(GN_inputs, emit);


%% Use the tblut retrieval as the initial guess for the hyperspectral retrieval

tic
% --------------------------------------------------------------
% ---------------- Retrieve Vertical Profile! ------------------
% --------------------------------------------------------------

[GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_4EMIT_top_bottom_ver2(GN_inputs, emit, spec_response, folder_paths);

disp([newline, 'Multispectral retrieval took ', num2str(toc), 'seconds to run', newline])

% --------------------------------------------------------------
% --------------------------------------------------------------


%%
% ----------------------------------------------
% ------------ SAVE OUTPUT STRUCTURE -----------
% ----------------------------------------------

% Save the version without an measurement uncertainty. Then we can add
% uncertainty and save the new file



if exist(folder_paths.saveOutput_filename, 'file')==2
    % append
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "pixels2use", "folder_paths", '-append');

else
    save(folder_paths.saveOutput_filename, "GN_outputs", "GN_inputs", "pixels2use", "folder_paths");

end


%% Make plot of the retrieved profile

plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)


