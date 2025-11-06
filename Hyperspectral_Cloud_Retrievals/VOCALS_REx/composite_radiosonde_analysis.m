%% Analyze composite soundings from VOCALS-REx


% By Andrew John Buggee
%%

clear variables

% Read all the file names

which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true


    folderpath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/allSoundings_composite_5mb/'];

    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

elseif strcmp(which_computer,'andrewbuggee')==true


    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/radiosonde/'];


    % ***** Define the ensemble profiles folder *****

    vocalsRexFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];


end



% ---------------------
% define the ensemble filename
ensemble_profiles = 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat';


%% Load the ensemble profiles

load([vocalsRexFolder, ensemble_profiles]);


%% What do I need to preform this analysis? 

% Ideally I would have a radiosonde measurement close to each vertical
% droplet profile, thus I could estimate the water vapor concentration
% profile for each droplet profile measured and therefore estimate
% statistics of the above-cloud precipitable water amount. 

% To do this, let's step through each ensemble profile, find the radiosonde
% closest to the sampled cloud, and store the radiosonde sampled vertical
% profiles

%% Step through each profile and find the radiosonde closest in space and time

for nn = 1:length(ensemble_profiles)

    % the composite radiosonde data is saved with the day it was recorded
    % on. So first find the radiosonde file for the day the current
    % ensemble profile was recorded

    % Extract the date from the ensemble profile filename
    dateStr = extractBetween(ensemble_profiles, 'from_', '_LWC-threshold');
    dateStr = char(dateStr); % Convert to char array for file naming

    % Construct the corresponding radiosonde filename
    radiosonde_filename = sprintf('radiosonde_%s.mat', dateStr);
    
    % Load the radiosonde data
    load(fullfile(folderpath, radiosonde_filename));






end

