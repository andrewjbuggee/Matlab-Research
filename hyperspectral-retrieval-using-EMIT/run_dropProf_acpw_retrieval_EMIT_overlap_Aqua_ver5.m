%% Run EMIT retireval SLURM function


function [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] =...
    run_dropProf_acpw_retrieval_EMIT_overlap_Aqua_ver5(folder_paths, print_status_updates, print_libRadtran_err,...
    plot_figures, save_figures)


which_computer = folder_paths.which_computer;


%% Define the folder of the coincident data set between EMIT and Aqau

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------


% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    % add folders to the path
    addpath(genpath('/projects/anbu8374/Matlab-Research'));
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/'));
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/'));


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

if plot_figures == true

    figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
    hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
    hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')
    hold on; geoscatter(amsr.geo.Latitude(:), amsr.geo.Longitude(:), 10, 'k.')

end

%% Plot the pixel footprints on the Earth to see the overlap
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
    % [fig1, ax1] = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);
    [fig1, ax1] = plot_instrument_footprints_4(modis, emit, airs, amsr, overlap_pixels, options);

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

        warning([newline, 'EMIT pixel ', num2str(pp), ' is masked out. Moving to pixel ', num2str(pp+1), '...', newline])

        continue

    else

        if pp>1

            clear GN_inputs GN_outputs tblut_retrieval acpw_retrieval

        end

        disp([newline, 'Retrieving Profile for pixel ', num2str(pp), '...', newline])

    end



    % *** Retrieve r_top, r_bot, tau_c, and acpw ***
    % *** CURRENT FORWARD MODEL UNCERTAINTIES CONSIDERED ***
    % (1) Adiabatic droplet profile assumption
    % (2) Cloud top height assumption
    
    % [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = retrieve_dropProf_acpw_EMIT_Aqua_singlePix_ver1(emit,...
    %     modis, airs, overlap_pixels,...
    %     folder_paths, print_libRadtran_err, print_status_updates, pp);



    % *** Retrieve r_top, r_bot, tau_c, and acpw ***
    % *** CURRENT FORWARD MODEL UNCERTAINTIES CONSIDERED ***
    % (1) Adiabatic droplet profile assumption
    % (2) Cloud top height assumption
    % (3) Droplet distribution effective variance assumption

    [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = retrieve_dropProf_acpw_EMIT_Aqua_singlePix_ver2(emit,...
        modis, airs, overlap_pixels,...
        folder_paths, print_libRadtran_err, print_status_updates, pp);







end









end
