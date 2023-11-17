% ----- RUN and ANALYZE MODIS L1N Data ----


clear variables;
scriptPlotting_blk;
% By Andrew J. Buggee
%% ---- Read Data ----

% Define data folders and files that you'd like to read

folderName = './MODIS_data/2021_08_25/L1B_calibrated_radiances/';
fileName = 'MOD02HKM.A2021237.1920.006.2021238072711.hdf';

geoFolderName = './MODIS_data/2021_08_25/GeoLocation_data/';
geoFileName = 'MOD03.A2021237.1920.006.2021238012056.hdf';


% extract data from hdf files
[EV] = readMODIS_L1B_data(folderName,fileName);
[sensor,solar,geo] = readMODIS_geolocation(geoFolderName,geoFileName);



%% -- plot bands --

plotRes = '250'; % which resolution to plot - '250', '500', or '1000'
quantity = 'reflectance'; % which quantity to plot - 'radiance' or 'reflectance'

% inputs for the plot function are: (data,
% sensor_parameters,location_params,plot_resolution,plot_quantity)

plotL1B_data(EV,geo,plotRes,quantity)


