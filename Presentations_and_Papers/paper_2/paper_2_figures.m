%% Figures for Paper 2

% By Andrew John Buggee

%%  Plot droplet profile retrievals using simualted HySICS measurements with vocals-Rex in-situ measurements





% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



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



plt_idx = 73;

plot_retrieved_prof_with_inSitu_paper2(folder_paths.retrieval, filenames_retrieval(plt_idx).name)


%% Plot the new figure with cloud top height unceratinty included in the covariance matrix
 
clear variables


folder_path = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    'HySICS/Droplet_profile_retrievals/paper2_variableSweep/test_logSpace_newCov_with_VR_inSitu_meas/'];

filename = ['dropletRetrieval_HySICS_636bands_0.3%_uncert_vocalsRex_recorded_',...
    '31-Oct-2008_9.0864UTC_vza_4_vaz_257_sza_31_saz_96_sim-ran-on-23-Dec-2025_1.mat'];

plot_retrieved_prof_with_inSitu_paper2(folder_path, filename)
