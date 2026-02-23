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
    retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_6';


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

% ----------------------------------------
% -- For Retrieval Results from Take 4 ---
% ----------------------------------------
% profile_indexes for paper = [1, 4, ]
plt_idx = 8;
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
