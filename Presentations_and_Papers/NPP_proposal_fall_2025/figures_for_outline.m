%% Plots for NPP proposal

% By Andrew John Buggee

%% Create plot of droplet profiles from VOCALS-REx that show an increase in cloud droplet effective radius at cloud top

clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true


    % Location of ensemble data
    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_10-Nov-2025'])



elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Location of ensemble data
    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
        'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat'])

    % --- all profiles, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    %     'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_02-Nov-2025.mat']);



end



% idx_2Plot = [2 ,6, 7, 9, 10, 11, 15, 16, 19, 20, 21, 28, 37, 43, 52, 57];
idx_2Plot = [2 ,6, 7, 9, 10, 11, ];



n_subPlots = 6;

ax = [];

for nn = 1:n_subPlots

     plot_LWC_and_re_and_Nc_vs_altitude(ensemble_profiles, idx_2Plot(nn), false)


end