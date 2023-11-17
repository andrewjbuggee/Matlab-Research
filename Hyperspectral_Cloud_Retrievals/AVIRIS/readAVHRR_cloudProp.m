%% Read AVHRR Cloud Properties Data

% This code will read AVHRR netCDF files, retreive clour propeties, sptial
% and temporal information associated with it, and output a struture. 

% By Andrew John Buggee

%% --- SETUP ---

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true
    
    folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/AVIRIS_Classic/AVHRR_Data/';
    fileName = 'patmosx_v05r03_NOAA-19_des_d20150205_c20170603.nc';
    
elseif strcmp(computer_name,'andrewbuggee')==true
    
    error('No Folders or Files yet!')   
end

%% --- READ IN THE DATA ---


% read in effective radius values
re = ncread([folderName,fileName], 'cld_reff_dcomp');
re_unc = ncread([folderName,fileName], 'cld_reff_dcomp_unc');

% read in effective radius values
tau = ncread([folderName,fileName], 'cld_opd_dcomp');
tau_unc = ncread([folderName,fileName], 'cld_opd_dcomp_unc');


% read in the time at which the data was taken
t = ncread([folderName,fileName], 'time');

% read in the lat and long of the data
lat = ncread([folderName,fileName], 'latitude');
long = ncread([folderName,fileName], 'longitude');

% create mesh grid
[Lat, Long] = meshgrid(lat,long);

%% --- Plot the Data ---

% Let's view the data on a geoplot

figure; geoscatter(Lat(:), Long(:));




