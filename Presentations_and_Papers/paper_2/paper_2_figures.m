%% Figures for Paper 2

% By Andrew John Buggee

%% Plot the HySICS retrieval along with the in-situ measurement

% ** only considering re_profile uncertainty **

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if strcmp(filenames_retrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];


% ------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_3 ---
% ------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
%plt_idx = 17;
% ------------------------------------------------------------


% -------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1 ---
% -------------------------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
plt_idx = 4;
% ------------------------------------------------------------


fig1 = plot_retrieved_prof_with_inSitu_paper2(folder_paths.retrieval, filenames_retrieval(plt_idx).name);



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
    num2str(plt_idx), '.fig']);


% save .png with 400 DPI resolution
% remove title
title('');
exportgraphics(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
    num2str(plt_idx), '.jpg'],'Resolution', 400);
% -------------------------------------
% -------------------------------------




%% Plot true color image from EMIT with overlapping footprints from the Aqua instruments


clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders(2);
which_computer = folder_paths.which_computer;


% Would you like to print status updates and/or the libRadtran error file?


plot_figures = true;

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
    % folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

    folder_paths.coincident_dataFolder = '2023_9_16_T191130_1/';

end








% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;   % cloud fraction
criteria.cld_tau_min = 3;   % cloud optical depth
criteria.cld_tau_max = 30;   % cloud optical depth
criteria.H = 0.1;         % horizontal inhomogeneity index





% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - Then compute: overlap.time_difference_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))

[overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
    criteria, plot_figures);

% ** If there aren't any pixels found ... **
% Increase the horizontal inhomogeneity index

while isempty(overlap_pixels.modis.linear_idx) == true

    disp([newline, 'No overlaping pixels that meet defined criteria. Increasing H index....', newline])
    criteria.H = criteria.H + 0.25;         % horizontal inhomogeneity index

    % recompute
    [overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
        criteria, plot_figures);

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

    % Values for 2023_9_16_T191130_1 scene
    % options.latlim = [-25.7, -25.46];  % Only show -30° to -20° latitude
    % options.lonlim = [-71.15, -70.8];  % Only show -75° to -65° longitude

    % ** Plot with RGB Image **
    % fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
    % fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    % [fig1, ax1] = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);
    [fig1, ax1] = plot_instrument_footprints_4(modis, emit, airs, amsr, overlap_pixels, options);
    % [fig1a, ax1a] = plot_instrument_footprints_4(modis, emit, [], amsr, overlap_pixels, options);


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


        % save .png with 500 DPI resolution
        % remove title
        ax1.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene with MODIS context - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

        f = gcf;
        ax1a.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene with MODIS context - Zoom In version - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

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



%% 

% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);



%% Plot EMIT retrievals!

% clear variables





% Load simulated measurements
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(whatComputer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    data_path.stored_retrievals = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/Droplet_profile_retrievals/',...
        'Paper_2/Droplet_profile_retrievals_take1/'];



    % ------------------------------------------
    % **********  2023_9_16_T191130  ***********
    % ------------------------------------------
    % 2023_9_16_T191130 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_1_ran-on-08-Jan-2026_rev1.mat';
    
    % 2023_9_16_T191130 pixel 2
    data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_2_ran-on-09-Jan-2026_rev1.mat';



    % ------------------------------------------
    % **********  2023_9_16_T191118  ***********
    % ------------------------------------------
    % 2023_9_16_T191118 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_1_ran-on-08-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 2 -
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_2_ran-on-09-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 3 - 
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_3_ran-on-09-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 4 - 
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_4_ran-on-09-Jan-2026_rev1.mat';


end

% make sure the MODIS, AIRS and AMSR-E data match the EMIT data used for
% the above retrieval
if strcmp(folder_paths.coincident_dataFolder(1:end-3),...
        extractBetween(data_path.retrieval_name, "Retrieval_", "_pixel")) == false

    error([newline, 'The MODIS, AIRS and AMSR-E data dont match the EMIT data time!', newline])

end

load([data_path.stored_retrievals, data_path.retrieval_name])

% what pixel number is this?
pixel_num_2Plot = str2double(extractBetween(data_path.retrieval_name, "pixel_", "_ran"));

% Make plot of the retrieved profile

% plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
fig3 = plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR_perPixel(GN_outputs, GN_inputs, modis, [], amsr, pixel_num_2Plot);

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
% % save .png with 400 DPI resolution
% % remove title
% exportgraphics(fig3,[folderpath_figs,...
%     'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);



%%