%% ----- READ IN MODIS DATA -----

% this function will read in L1A MODIS data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf


% ---- Description of the different fields within the HDF File -----




% By Andrew J. Buggee
%%

function [EV,geo,solar,sensor] = readMODIS_L1A_data(folderName,fileName,geoFileName)
%% ---- Read in File Info -----




%data_info = hdfinfo([folderName,fileName]);

%% --- Read in a data set at a specifc resolution ---

% for now we only read in Earth-viewing data at two of 3 resolutions: 250m
% and 500m. The 1km resolution data may be used later on. But this consists
% of 29 spectral bands. So to save loacl memory we wont load these yet



    
    earthView250 = hdfread([folderName,fileName],'EV_250m'); % 250 meter resolution data which carries spectral bands 1 and 2

    earthView500 = hdfread([folderName,fileName],'EV_500m'); % 500 meter resolution data which carries spectral bands 3,4,5,6 and 7

    

%     earthView1km_day = hdfread([folderName,fileName],'EV_1km_day');
%     earthView1km_night = hdfread([folderName,fileName],'EV_1km_night');
    
    % load the geolocation data from the MOD03 geolocation hdf file
    % these are the lat-long positions of MODIS pixels on Earths surface
    geo.lat = hdfread([folderName,geoFileName],'Latitude');
    geo.long = hdfread([folderName,geoFileName],'Longitude');
    

    % load solar position data
    solar.azimuth = hdfread([folderName,geoFileName],'SolarAzimuth');
    solar.zenith = hdfread([folderName,geoFileName],'SolarZenith');
    
    % load satellite position data
    sensor.height = hdfread([folderName,geoFileName],'Height');
    sensor.range = hdfread([folderName,geoFileName],'Range');
    sensor.azimuth = hdfread([folderName,geoFileName],'SensorAzimuth');
    sensor.zenith = hdfread([folderName,geoFileName],'SensorZenith');
    

% split the 250 meter earth view data in (row,column,band) format

numRows = size(earthView250,1); % across track direction, taken using scan mirror
numCols = size(earthView250,3); % frames?
numBands = size(earthView250,2); % bands

EV250 = zeros(numRows,numCols,numBands);

for ii = numBands
    
    EV250(:,:,ii) = reshape(earthView250(:,ii,:),numRows,numCols);
    
end

EV.m250 = EV250;

% split the 500 meter earth view data in (row,column,band) format

numRows = size(earthView500,1); % across track direction, taken using scan mirror
numCols = size(earthView500,3); % frames?
numBands = size(earthView500,2); % bands

EV500 = zeros(numRows,numCols,numBands);

for ii = numBands
    
    EV500(:,:,ii) = reshape(earthView500(:,ii,:),numRows,numCols);
    
end

EV.m500 = EV500;



%% ---- MODIS Parameters -----

sensor.modisBands = [645,858.5, 469, 550, 1240, 1640, 2130];

end





