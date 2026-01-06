%% Optimal estimation of a vertical droplet profile using EMIT data

% Retrieve a droplet profile from EMIT data
% Use spectral reflectances outside of water vapor bands


% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


%% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders(2);
which_computer = folder_paths.which_computer;


%% Would you like to print status updates and/or the libRadtran error file?

print_status_updates = true;

print_libRadtran_err = false;


%% Define the folder of the coincident data set between EMIT and Aqau

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------


% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folder_paths.coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    folder_paths.coincident_dataFolder = '2024-09-12/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    % folder_paths.coincident_dataFolder = '2024-09-12/';

    % folder_paths.coincident_dataFolder = '2024_05_17-T1835/';

    % EMIT pixels masked out
    % folder_paths.coincident_dataFolder = '2023_9_16_T191106_2/';

    % 11 Pixels with H less than 1.6
    % folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

    % 2 Pixels with H less than 1.1
    % folder_paths.coincident_dataFolder = '2023_9_16_T191130_1/';

    % 14 Pixels with H less than 1.35
    % folder_paths.coincident_dataFolder = '2023_9_16_T191142_1/';

    % 2 Pixels with H less than 1.35
    folder_paths.coincident_dataFolder = '2024_1_13_T194658_1/';

elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));
    addLibRadTran_paths;

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];


    % 2 Pixels with H less than 1.35
    folder_paths.coincident_dataFolder = '2024_1_13_T194658_1/';

end






%% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;   % cloud fraction
criteria.cld_tau_min = 3;   % cloud optical depth
criteria.cld_tau_max = 30;   % cloud optical depth
criteria.H = 0.1;         % horizontal inhomogeneity index

% plot flag
plot_data = false;

% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - Then compute: overlap.time_difference_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))

[overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
    criteria, plot_data);

% ** If there aren't any pixels found ... **
% Increase the horizontal inhomogeneity index

while isempty(overlap_pixels.modis.linear_idx) == true

    disp([newline, 'No overlaping pixels that meet defined criteria. Increasing H index....', newline])
    criteria.H = criteria.H + 0.25;         % horizontal inhomogeneity index

    % recompute
    [overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
        criteria, plot_data);

    if isempty(overlap_pixels.modis.linear_idx) == false

        % print the H value used
        disp([newline, 'Horizontal Inhomogeneity Index - H = ', num2str(criteria.H), newline])

    end

end

%% Plot all three swaths

figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')
hold on; geoscatter(amsr.geo.Latitude(:), amsr.geo.Longitude(:), 10, 'k.')

%% Plot the pixel footprints on the Earth to see the overlap
% Add an RGB true color image for context
clear options
options.use_radiance = false;
[rgb_img, rgb_lat, rgb_lon] = create_modis_true_color(modis, options);

options.show_rgb = true;
options.rgb_image = rgb_img;
options.rgb_lat = rgb_lat;
options.rgb_lon = rgb_lon;
options.latlim = [-37, -27];  % Only show -30° to -20° latitude
options.lonlim = [-80, -70];  % Only show -75° to -65° longitude

% ** Plot with RGB Image **
% fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);

% ** Plot without RGB Image **
options.show_rgb = false;
fig2 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);


%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);




%% Set libRadtran INP directory

folder_paths.libRadtran_inp = [folder_paths.libRadtran_inp, 'EMIT_', folder_paths.coincident_dataFolder(1:end-1), '_',...
    folder_paths.L1B_fileName_emit{1}(27:30),'/'];


% If the folder path doesn't exit, create a new directory
if ~exist(folder_paths.libRadtran_inp, 'dir')

    mkdir(folder_paths.libRadtran_inp)

end


%%   Delete old files?

% First, delete files in the HySICS folder
delete([folder_paths.libRadtran_inp, '*.INP'])
delete([folder_paths.libRadtran_inp, '*.OUT'])

% delete old wc files
delete([folder_paths.libRadtran_water_cloud_files, '*.DAT'])

% delete old atm files
delete([folder_paths.atm_folder_path, '*.DAT'])

% delete old mie files
delete([folder_paths.libRadtran_mie_folder, '*.INP'])
delete([folder_paths.libRadtran_mie_folder, '*.OUT'])






%% RUN FOR LOOP OVER ALL PIXELS

for pp = 1:length(overlap_pixels.modis.linear_idx)


    %% Check to see if all chosen EMIT pixels have been masked out

    if all(isnan(emit.radiance.measurements(:, pp)))

        warning([newline, 'EMIT pixels are masked out. Moving to pixel ', num2str(pp), '...', newline])

        continue

    else

        if pp>1

            clear GN_inputs GN_outputs tblut_retrieval acpw_retrieval

        end

        disp([newline, 'Retrieving Profile for pixel ', num2str(pp), '...', newline])

    end
    

    [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = run_retrieval_dropProf_acpw_EMIT_Aqua_singlePix_ver1(emit,...
            modis, airs, overlap_pixels,...
            folder_paths, print_libRadtran_err, print_status_updates, pp);





    %% Make plot of the retrieved profile

    % plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
    % plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR(GN_outputs, GN_inputs, modis, [], amsr)



end


