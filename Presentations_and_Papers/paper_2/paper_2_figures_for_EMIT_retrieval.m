%% Figures for Paper 2

% By Andrew John Buggee



%% Plot true color image from EMIT with overlapping footprints from the Aqua instruments


clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders(2);
which_computer = folder_paths.which_computer;


% Would you like to print status updates and/or the libRadtran error file?


plot_figures = false;

save_figures = false;


% Define the folder of the coincident data set between EMIT and Aqau

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
    folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';


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


    folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

end








% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

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




% Plot all three swaths

if plot_figures == true

    figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
    hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
    hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')
    hold on; geoscatter(amsr.geo.Latitude(:), amsr.geo.Longitude(:), 10, 'k.')

end




% Plot the pixel footprints on the Earth to see the overlap
% Add an RGB true color image for context

if plot_figures == true

    clear options
    % options.use_radiance = false;
    % options.rgb_image_type = 'modis';
    % [rgb_img, rgb_lat, rgb_lon] = create_modis_true_color(modis, options);

    options.convert_to_reflectance = false;
    options.rgb_image_type = 'emit';
    [rgb_img, rgb_lat, rgb_lon, band_indices] = create_emit_true_color(emit, options);


    options.show_rgb = true;
    options.rgb_image = rgb_img;
    options.rgb_lat = rgb_lat;
    options.rgb_lon = rgb_lon;
    % options.latlim = [-30, -20];  % Only show -30° to -20° latitude
    % options.lonlim = [-80, -67];  % Only show -75° to -65° longitude

    % ** Plot with RGB Image **
    % fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
    % fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    [fig1, ax1] = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);

    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    % save .fig file
    if save_figures==true

        if strcmp(which_computer,'anbu8374')==true
            error(['Where do I save the figure?'])
        elseif strcmp(which_computer,'andrewbuggee')==true
            folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
        end
        f = gcf;
        saveas(f,[folderpath_figs,'EMIT Scene with MODIS context - ', folder_paths.coincident_dataFolder(1:end-1), '.fig']);


        % save .png with 400 DPI resolution
        % remove title
        ax1.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene with MODIS context - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 400);

    end
    % -------------------------------------
    % -------------------------------------


    % ** Plot without RGB Image **
    options.show_rgb = false;
    % fig2 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    fig2 = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);

    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    if save_figures==true

        % save .fig file
        if strcmp(which_computer,'anbu8374')==true
            error(['Where do I save the figure?'])
        elseif strcmp(which_computer,'andrewbuggee')==true
            folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
        end
        % remove title
        ax1.Title.String = '';
        f = gcf;
        saveas(f,[folderpath_figs,'EMIT Scene and Aqua instrument overlap without MODIS context - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.fig']);


        % save .png with 400 DPI resolution
        % remove title
        ax1.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene and Aqau instrument overlap without MODIS context - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 400);

    end
    % -------------------------------------
    % -------------------------------------

end


%% Plot VOCALS-REx CDP and radiosonde data

% plot VOCALS-REx in-situ measured re, lwc, and Nc as a function of
% altitude, along with the radiosonde measured temperature and RH


clear variables


% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    % ***** Define the radiosonde profiles folder *****
    folderpath_radioSonde = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/paper2_prior_stats/'];
    
    radiosonde_profs = 'radiosonde_profiles_for_paper2_measurements_closest_radiosonde_in_time_27-Jan-2026.mat';


    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    % define the ensemble filename
    % profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025.mat';
    profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat';

elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath_radioSonde = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];


    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end


% Load the ensemble profiles and the radisonde profiles

airborne = load([vocalsRexFolder, profiles]);
radSonde = load([folderpath_radioSonde, radiosonde_profs]);


% define the vert prof to plot: [1, 69]

plt_idx = 3;

plot_VOCALS_insitu_re_lwc_nc_and_radioSonde(airborne.ensemble_profiles, radSonde, plt_idx, false)




%% Plots EMIT retrieval results

clear variables

% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    retrieval_directory = '/Users/anbu8374/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4/';

    coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------
    % -------------------------------------

    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4/';
    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_5';
    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_6';

    retrieval_directory = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7'];


    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';



end



% Grab filenames in drive
filenames_retrieval = dir(retrieval_directory);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == false

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];

% ----------------------------------------
% -- For Retrieval Results from Take 4 ---
% ----------------------------------------
% profile_indexes for paper = [5, 12, ]
plt_idx = 36;
% ------------------------------------------------------------

load([filenames_retrieval(plt_idx).folder, '/', filenames_retrieval(plt_idx).name])




% ----------------------------------------
% *** Extract the pixel number ***
% ----------------------------------------
pixel_num = str2double(extractBetween([filenames_retrieval(plt_idx).folder, '/', filenames_retrieval(plt_idx).name],...
    'pixel_', '_'));



% ----------------------------------------
% *** Load MODIS, AIRS and AMSR-E data ***
% ----------------------------------------

% Load EMIT data
[emit, folder_paths.L1B_fileName_emit] = retrieveEMIT_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load Aqua/MODIS Data
[modis, ~] = retrieveMODIS_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load AIRS data
airs = readAIRS_L2_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load AMSR-E/2 data
amsr = readAMSR_L2_data([coincident_dataPath, folder_paths.coincident_dataFolder]);
% ----------------------------------------




% ----------------------------------------
% Remove data that is not needed
% ----------------------------------------

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);




% ----------------------------------------------------
% Compute the above cloud precipitable water from AIRS
% ----------------------------------------------------
if isfield(airs, 'acpw') == false

    assumed_cloudTopHeight = GN_inputs.RT.z_topBottom(1)*1e3;    % meters - cloud top height used by libRadtran

    % Compute the above cloud precipitable water from AIRS data
    airs = convert_AIRS_prof_2_mass_density(airs, atm_data_directory,...
        pixel_num, overlap_pixels, [], false, assumed_cloudTopHeight);

end




% -------------------------------------------------------
% -------------------------------------------------------
use_new_LWP_calc = true;

if use_new_LWP_calc == true

% ** Compute new updated LWP calc ***
re_profile = create_droplet_profile2([GN_outputs.retrieval(1,end), GN_outputs.retrieval(2,end)],...
    GN_inputs.RT.z, 'altitude', GN_inputs.model.profile.type);


% define the z vector
z = linspace(GN_inputs.RT.z_topBottom(2), GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

% define the z midpoint at each layer and normalize it!
z_norm = z - z(1);
z_norm_mid = (diff(z_norm)/2 + z_norm(1:end-1));


% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
ext_bulk_coeff_per_LWC = zeros(length(re_profile), 1);

for rr = 1:length(re_profile)

    mie_radius = [re_profile(rr), re_profile(rr), 0];    % microns

    size_distribution = {'gamma', GN_inputs.RT.distribution_var(rr)};           % droplet distribution

    % Create a mie file
    [input_filename, output_filename] = write_mie_file('MIEV0', 'water',...
        mie_radius, 500, size_distribution, 'verbose', rr, round(re_profile(rr), 4), mie_folder_path);

    % run the mie file
    [~] = runMIE(mie_folder_path, input_filename,output_filename, which_computer);

    % Read the output of the mie file
    [mie,~,~] = readMIE(mie_folder_path, output_filename);

    ext_bulk_coeff_per_LWC(rr) = mie.Qext;       % km^-1 / (cm^3 / m^3)

end


% ** Assuming liquid water content increases linearly with depth **
z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
dz_km = z(2) - z(1);           % kilometers

slope = GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

% solve for the linear liquid water content profile
%lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


lwp_newCalc = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2

GN_outputs.LWP_newCalc = lwp_newCalc;

end

% -------------------------------------------------------
% -------------------------------------------------------





% plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
fig3 = plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR_perPixel(GN_outputs, GN_inputs, modis,...
    airs, amsr, pixel_num, overlap_pixels, use_new_LWP_calc);

% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(which_computer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(which_computer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig3,[folderpath_figs,'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1), '.fig']);


% save .png with 500 DPI resolution
% remove title
exportgraphics(fig3,[folderpath_figs,...
    'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1),...
    'pixelNum-', num2str(pixel_num), '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------





%% Plots the mean EMIT retrieval over N pixels within a single MODIS pixel

clear variables

% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    retrieval_directory = '/Users/anbu8374/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4/';

    coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------
    % -------------------------------------

    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4/';
    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_5';
    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_6';

    retrieval_directory = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7'];


    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';



end



% Grab filenames in drive
filenames_retrieval = dir(retrieval_directory);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == false

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];

retrieval_output = cell(7, 1);   % There are seven directories 

date_str = cell(length(filenames_retrieval), 1);

% Step through the file names and store the retirevals over different
% pixels for the same scence
for nn = 1:length(filenames_retrieval)

    date_str(nn) = extractBetween(filenames_retrieval(nn).name, 'dropRetrieval_', '_EMIT');
end

% Find the number of unique date strings
date_str = unique(date_str);

% Now step through each filename again and store data for each unique day
for nn = 1:length(filenames_retrieval)

    date_str_2save = extractBetween(filenames_retrieval(nn).name, 'dropRetrieval_', '_EMIT');

    for mm = 1:length(date_str)

        idx_2save = [];
        if strcmp(date_str{mm}, date_str_2save) == true

            idx_2save = mm;

            break   

        end

    end

    % Store the data
    load([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name])
    retrieval_output{idx_2save} = [retrieval_output{idx_2save}, GN_outputs.retrieval(:, end)];



end


        



% ----------------------------------------
% -- For Retrieval Results from Take 4 ---
% ----------------------------------------
% profile_indexes for paper = [1, 4, ]
plt_idx = 1;
% ------------------------------------------------------------

load([filenames_retrieval(plt_idx).folder, '/', filenames_retrieval(plt_idx).name])



% ----------------------------------------
% *** Extract the pixel number ***
% ----------------------------------------
pixel_num = str2double(extractBetween([filenames_retrieval(plt_idx).folder, '/', filenames_retrieval(plt_idx).name],...
    'pixel_', '_'));




% ----------------------------------------
% *** Load MODIS, AIRS and AMSR-E data ***
% ----------------------------------------

% Load EMIT data
[emit, folder_paths.L1B_fileName_emit] = retrieveEMIT_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load Aqua/MODIS Data
[modis, ~] = retrieveMODIS_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load AIRS data
airs = readAIRS_L2_data([coincident_dataPath, folder_paths.coincident_dataFolder]);

% Load AMSR-E/2 data
amsr = readAMSR_L2_data([coincident_dataPath, folder_paths.coincident_dataFolder]);
% ----------------------------------------




% ----------------------------------------
% Remove data that is not needed
% ----------------------------------------

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);




% ----------------------------------------------------
% Compute the above cloud precipitable water from AIRS
% ----------------------------------------------------
if isfield(airs, 'acpw') == false

    assumed_cloudTopHeight = GN_inputs.RT.z_topBottom(1)*1e3;    % meters - cloud top height used by libRadtran

    % Compute the above cloud precipitable water from AIRS data
    airs = convert_AIRS_prof_2_mass_density(airs, atm_data_directory,...
        pixel_num, overlap_pixels, [], false, assumed_cloudTopHeight);

end


% plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
fig3 = plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR_perPixel(GN_outputs, GN_inputs, modis,...
    airs, amsr, pixel_num, overlap_pixels);

% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(which_computer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(which_computer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
% end
% saveas(fig3,[folderpath_figs,'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1), '.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% exportgraphics(fig3,[folderpath_figs,...
%     'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------


%% Look at all of the EMIT retrievals and see how many advance beyone the initial guess


%% Plots EMIT retrieval results


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    retrieval_directory = '/Users/anbu8374/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4/';

    coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------
    % -------------------------------------
    
    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/Droplet_profile_retrievals/Paper_2/take_4';

    % retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_5';

    % retrieval_directory = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/',...
    %     'Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/Droplet_profile_retrievals/take_7'];

    retrieval_directory = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7'];


    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';



end



% Grab filenames in drive
filenames_retrieval = dir(retrieval_directory);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == true

        load([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name])

        disp([newline, 'File: ', num2str(nn), newline])
        disp(['Filename: ', filenames_retrieval(nn).name, newline])
        disp(['     Retrieval Iterations: ', newline])
        disp(['                   re_top: ', num2str(GN_outputs.retrieval(1, :)), newline])
        disp(['                   re_bot: ', num2str(GN_outputs.retrieval(2, :)), newline])
        disp(['                    Tau_c: ', num2str(GN_outputs.retrieval(3, :)), newline])
        disp(['                     acpw: ', num2str(GN_outputs.retrieval(4, :)), newline])
        disp(['----------------------------------------', newline])
        disp(['       RSS convergence limit: ', num2str(GN_inputs.convergence_limit), newline])
        disp(['RSS using final state vector: ', num2str(GN_outputs.rss_residual(end)), newline])
        
  

    end

end



