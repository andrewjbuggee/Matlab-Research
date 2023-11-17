%% Find the intersection between Vocals-Rex flight path and the MODIS pixels


% By Andrew John Buggee

%% Load MODIS and VOCALS-REx data



% Load modis data and create input structure
clear variables

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';

    % ----------------------------------------
    % ***** Define the VOCALS-REx Folder *****
    % ----------------------------------------

    % ----- November 9th data -----
    vocalsRexFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/2008_11_09/';
    vocalsRexFile = 'RF11.20081109.125700_213600.PNI.nc';



    % ----- November 11 data -----
%     vocalsRexFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/2008_11_11/';
%     vocalsRexFile = 'RF12.20081111.125000_214500.PNI.nc';



elseif strcmp(computer_name,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';




    % ----------------------------------------
    % ***** Define the VOCALS-REx Folder *****
    % ----------------------------------------



    % ----- November 9 data -----
    vocalsRexFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/2008_11_09/';
    vocalsRexFile = 'RF11.20081109.125700_213600.PNI.nc';


end




[modis,L1B_fileName] = retrieveMODIS_data(modisFolder);


% Load VOCALS-REx data
vocalsRex = readVocalsRex([vocalsRexFolder,vocalsRexFile]);

%% Find the vocals-rex vertical profile closest in time to the MODIS data 
% only keep the vocals rex data closest to MODIS

% ----- Find all vertical profiles within VOCALS-REx data ------
% find all vertical profiles in the vocals-REx data set. Vertical profiles
% qualify if the total number concentration is above 10^0 and relatively
% flat, the plane is climbing in altitude, and clears the cloud deck from
% bottom to top. If this is all true, save the vocals rex data
lwc_threshold = 0.03;           % g/m^3
stop_at_max_lwc = true;         % truncate profile after the maximum lwc value

% Lets only keep the vocalsRex data we need
% Time is measured in seconds since the startTime

vocalsRex = cropVocalsRex2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, modis);

%%  Plot MODIS pixels with real projecttion near VOCALS-REx profile

% The MODIS lat/long points represent the center of the 1km pixels
% Use the satellite range to find the boundaries of the pixels
% The swath seen on the ground is different for nadir pixels than it is for
% pixels at the endge of the focal plane array. That's because MODIS uses a
% scan mirror to obtain the full swath of measurements. THe pixels at the
% edge are sampling a larger section of the ground because they travel
% through a longer path length. How do I visualize the volumes that are
% sampled, and the projection of the pixels on the ground?

