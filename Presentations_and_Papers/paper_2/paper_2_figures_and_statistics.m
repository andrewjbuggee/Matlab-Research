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
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];

    % define the folder where retrievals are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];



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
% plt_idx = 4;
% ------------------------------------------------------------


% ---------------------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3 ---
% ---------------------------------------------------------------------------------------------
% profile_indexes for paper = [1,2,4, 6, 7,]
plt_idx = 2;
% ------------------------------------------------------------


fig1 = plot_retrieved_prof_with_inSitu_paper2(folder_paths.retrieval, filenames_retrieval(plt_idx).name);



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% % save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
% end
% saveas(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
%     num2str(plt_idx), '.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
%     num2str(plt_idx), '.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------



%% Compute HySICS retrieval statistics



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

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    % % define the folder where retrievals are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_3/'];


    % define the folder where retrievals are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];

    addpath(folder_paths.retrieval);



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



con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3



% store the LWP retrieval
% store the TBLUT LWP estimate
% store the TBLUT LWP estimate with Wood-Hartmann adjustement
% store the in-situ LWP measurement
lwp_retrieval = zeros(size(filenames_retrieval));
lwp_tblut = zeros(size(filenames_retrieval));
lwp_tblut_WH = zeros(size(filenames_retrieval));
lwp_inSitu = zeros(size(filenames_retrieval));

lwp_newCalc = zeros(size(filenames_retrieval));


% store the ACPW retrieval
% store the true ACPW used in the forward model - measured by....
acpw_retrieval = zeros(size(filenames_retrieval));
acpw_inSitu = zeros(size(filenames_retrieval));


% store the optical depth retrieval
% store the in-situ measured optical depth
tauC_retrieval = zeros(size(filenames_retrieval));
tauC_inSitu = zeros(size(filenames_retrieval));



for nn = 1:length(filenames_retrieval)

    clear ds

    ds = load(filenames_retrieval(nn).name);

    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the liquid water paths
    lwp_retrieval(nn) = ds.GN_outputs.LWP;    % g/m^2

    % compute the LWP estimate using the TBLUT retrieval
    lwp_tblut(nn) = (2 * rho_h2o * (ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

    % ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
    lwp_tblut_WH(nn) = 5/9 * rho_h2o * ds.tblut_retrieval.minTau * (ds.tblut_retrieval.minRe/1e6); % g/m^2


    % What is the true LWP
    lwp_inSitu(nn) = ds.GN_inputs.measurement.lwp;   % g/m^2
    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % ** Compute new updated LWP calc ***
    re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
        ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


    % define the z vector
    z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

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

        size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

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

    %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
    slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

    % solve for the linear liquid water content profile
    %lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
    lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


    lwp_newCalc(nn) = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    acpw_retrieval(nn) = ds.GN_outputs.retrieval(end,end);    % mm

    % What is the true LWP
    acpw_inSitu(nn) = ds.GN_inputs.measurement.actpw;   % mm
    % -------------------------------------------------------
    % -------------------------------------------------------


    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    tauC_retrieval(nn) = ds.GN_outputs.retrieval(3,end);    %

    % What is the true LWP
    tauC_inSitu(nn) = ds.GN_inputs.measurement.tau_c;   %
    % -------------------------------------------------------
    % -------------------------------------------------------



end

% -------------------------------------------------------
% Compute statistics!!

% Let's compute the root-mean-square percent error
rms_err_lwp_hyperspectral = 100 * sqrt( mean( (1 - lwp_retrieval./lwp_inSitu).^2 ));  % percent
rms_err_lwp_tblut = 100 * sqrt( mean( (1 - lwp_tblut./lwp_inSitu).^2 ));  % percent
rms_err_lwp_tblut_WH = 100 * sqrt( mean( (1 - lwp_tblut_WH./lwp_inSitu).^2 ));  % percent

% using the new LWP calc!
rms_err_lwp_hyperspectral_newCalc = 100 * sqrt( mean( (1 - lwp_newCalc./lwp_inSitu).^2 ));  % percent


% Let's compute the average percent difference
avg_percent_diff_newCacl = mean( abs( 100 .* (1 - lwp_newCalc./lwp_inSitu) ));
avg_percent_diff_tblut = mean( abs( 100 .* (1 - lwp_tblut./lwp_inSitu) ));
avg_percent_diff_tblutWH = mean( abs( 100 .* (1 - lwp_tblut_WH./lwp_inSitu) ));







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
    folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

    % folder_paths.coincident_dataFolder = '2023_9_16_T191130_1/';

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



    % ** Plot with RGB Image **
    % fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
    % fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    % [fig1, ax1] = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);
    % [fig1, ax1] = plot_instrument_footprints_4(modis, emit, airs, amsr, overlap_pixels, options);


    % Values for 2023_9_16_T191130_1 scene
    % options.latlim = [-25.7, -25.46];  % Only show -30° to -20° latitude
    % options.lonlim = [-71.15, -70.8];  % Only show -75° to -65° longitude

    % Values for 2023_9_16_T191118_1 scene
    options.latlim = [-25.8, -25.5];  % Only show -30° to -20° latitude
    options.lonlim = [-71.53, -71];  % Only show -75° to -65° longitude

    [fig1a, ax1a] = plot_instrument_footprints_4(modis, emit, [], amsr, overlap_pixels, options);


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
            'EMIT Scene with MODIS context - Zoom In version - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

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
            folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

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

    addpath(data_path.stored_retrievals)



    % ------------------------------------------
    % **********  2023_9_16_T191130  ***********
    % ------------------------------------------
    % 2023_9_16_T191130 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_1_ran-on-08-Jan-2026_rev1.mat';

    % 2023_9_16_T191130 pixel 2
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_2_ran-on-09-Jan-2026_rev1.mat';



    % ------------------------------------------
    % **********  2023_9_16_T191118  ***********
    % ------------------------------------------
    % 2023_9_16_T191118 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_1_ran-on-08-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 2 -
    data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_2_ran-on-09-Jan-2026_rev1.mat';

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
if strcmp(which_computer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(which_computer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig3,[folderpath_figs,'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1),...
    '_pixel_num_', num2str(pixel_num_2Plot),'.fig']);


% save .png with 500 DPI resolution
% remove title
exportgraphics(fig3,[folderpath_figs,...
    'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1),...
    '_pixel_num_', num2str(pixel_num_2Plot), '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------







%% Compute EMIT-Aqua retrieval statistics



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

    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';

    % % define the folder where retrievals are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_3/'];


    % define the folder where retrievals are located
    % *** 2/13/2026 - Retrieval with overlapping EMIT/Aqua data
    %          
    % folder_paths.retrieval = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_5';


    % define the folder where retrievals are located
    % *** 2/15/2026 - Retrieval with overlapping EMIT/Aqua data             
    folder_paths.retrieval = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7';



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

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == false

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];



con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3



% store the LWP retrieval
% store the TBLUT LWP estimate
% store the TBLUT LWP estimate with Wood-Hartmann adjustement
% store the in-situ LWP measurement
lwp_retrieval = zeros(size(filenames_retrieval));
lwp_modis = zeros(size(filenames_retrieval));
lwp_modis_WH = zeros(size(filenames_retrieval));

lwp_newCalc = zeros(size(filenames_retrieval));


% store the ACPW retrieval
% store the true ACPW used in the forward model - measured by....
acpw_retrieval = zeros(size(filenames_retrieval));
acpw_modis = zeros(size(filenames_retrieval));
acpw_airs = zeros(size(filenames_retrieval));


% store the optical depth retrieval
% store the in-situ measured optical depth
tauC_retrieval = zeros(size(filenames_retrieval));
tauC_modis = zeros(size(filenames_retrieval));





for nn = 1:length(filenames_retrieval)

    clear ds emit modis airs

    ds = load([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name]);



    % ----------------------------------------
    % *** Extract the pixel number ***
    % ----------------------------------------
    pixel_num = str2double(extractBetween([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name],...
        'pixel_', '_'));



    % ----------------------------------------
    % *** Load MODIS, AIRS and AMSR-E data ***
    % ----------------------------------------

    % Load EMIT data
    % [emit, ~] = retrieveEMIT_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load Aqua/MODIS Data
    [modis, ~] = retrieveMODIS_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load AIRS data
    airs = readAIRS_L2_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load AMSR-E/2 data
    % amsr = readAMSR_L2_data([coincident_dataPath, folder_paths.coincident_dataFolder]);
    % ----------------------------------------

    % ----------------------------------------
    % Remove data that is not needed
    % ----------------------------------------
    % emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

    modis = remove_unwanted_modis_data(modis, ds.overlap_pixels.modis);

    airs = remove_unwanted_airs_data(airs, ds.overlap_pixels.airs);

    % amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);



    % Compute the above cloud precipitable water from AIRS data
    airs = convert_AIRS_prof_2_mass_density(airs, atm_data_directory,...
        pixel_num, ds.overlap_pixels, [], false, ds.GN_inputs.RT.z_topBottom(1)*1e3);



    % ** use the AIRS measurement closest to EMIT **
    unique_airs_pix = unique(ds.overlap_pixels.airs.linear_idx);
    unique_pix_idx_airs = zeros(1, length(ds.overlap_pixels.airs.linear_idx));
    for xx = 1:length(unique_pix_idx_airs)

        unique_pix_idx_airs(xx) = find(unique_airs_pix==ds.overlap_pixels.airs.linear_idx(xx));

    end


    % ** use the MODIS measurement closest to EMIT **
    unique_modis_pix = unique(ds.overlap_pixels.modis.linear_idx);
    unique_pix_idx_modis = zeros(1, length(ds.overlap_pixels.modis.linear_idx));
    for xx = 1:length(unique_pix_idx_modis)

        unique_pix_idx_modis(xx) = find(unique_modis_pix==ds.overlap_pixels.modis.linear_idx(xx));

    end



    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the liquid water paths
    lwp_retrieval(nn) = ds.GN_outputs.LWP;    % g/m^2

    % compute the LWP estimate using the TBLUT retrieval
    lwp_modis(nn) = (2 * rho_h2o *...
        (modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) )/1e6) *...
        modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) ) )/3; % g/m^2

    % ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
    lwp_modis_WH(nn) = 5/9 * rho_h2o *...
        (modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) )/1e6) *...
        modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) ); % g/m^2

    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % ** Compute new updated LWP calc ***
    re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
        ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


    % define the z vector
    z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

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

        size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

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

    %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
    slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

    % solve for the linear liquid water content profile
    %lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
    lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


    lwp_newCalc(nn) = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    acpw_retrieval(nn) = ds.GN_outputs.retrieval(end,end);    % mm

    % What is the true LWP
    acpw_modis(nn) = modis.vapor.col_nir( unique_pix_idx_modis(pixel_num) ) * 10;   % mm


    acpw_airs(nn) = airs.H2O.acpw_using_assumed_CTH( unique_pix_idx_airs(pixel_num) );
    % -------------------------------------------------------
    % -------------------------------------------------------


    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the retrieved optical depth
    tauC_retrieval(nn) = ds.GN_outputs.retrieval(3,end);    %

    % What is the MODIS optical depth
    tauC_modis(nn) = modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) );
    % -------------------------------------------------------
    % -------------------------------------------------------



end

% -------------------------------------------------------
% Compute statistics!!

% % Let's compute the root-mean-square percent error
% rms_err_lwp_hyperspectral = 100 * sqrt( mean( (1 - lwp_retrieval./lwp_inSitu).^2 ));  % percent
% rms_err_lwp_tblut = 100 * sqrt( mean( (1 - lwp_tblut./lwp_inSitu).^2 ));  % percent
% rms_err_lwp_tblut_WH = 100 * sqrt( mean( (1 - lwp_tblut_WH./lwp_inSitu).^2 ));  % percent
% 
% % using the new LWP calc!
% rms_err_lwp_hyperspectral_newCalc = 100 * sqrt( mean( (1 - lwp_newCalc./lwp_inSitu).^2 ));  % percent


% Let's compute the average percent difference for LWP
avg_percent_LWP_diff_newCacl_MODIS = mean( abs( 100 .* (1 - lwp_newCalc./lwp_modis) ));
avg_percent_LWP_diff_newCalc_MODIS_WH = mean( abs( 100 .* (1 - lwp_newCalc./lwp_modis_WH) ));

avg_percent_LWP_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - lwp_newCalc./lwp_modis) ));
avg_percent_LWP_diff_newCalc_MODIS_WH_noAbs = mean( ( 100 .* (1 - lwp_newCalc./lwp_modis_WH) ));



% Let's compute the average percent difference for ACPW
avg_percent_ACPW_diff_newCacl_MODIS = mean( abs( 100 .* (1 - acpw_retrieval./acpw_modis) ));
avg_percent_ACPW_diff_newCalc_AIRS = mean( abs( 100 .* (1 - acpw_retrieval./acpw_airs) ));

avg_percent_ACPW_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - acpw_retrieval./acpw_modis) ));
avg_percent_ACPW_diff_newCalc_AIRS_noAbs = mean( ( 100 .* (1 - acpw_retrieval./acpw_airs) ));

avg_percent_ACPW_diff_MODIS_AIRS_noAbs = mean( ( 100 .* (1 - acpw_airs./acpw_modis) ));


% Let's compute the average percent difference for optical thickness
avg_percent_tau_diff_newCacl_MODIS = mean( abs( 100 .* (1 - tauC_retrieval./tauC_modis) ));

avg_percent_tau_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - tauC_retrieval./tauC_modis) ));



