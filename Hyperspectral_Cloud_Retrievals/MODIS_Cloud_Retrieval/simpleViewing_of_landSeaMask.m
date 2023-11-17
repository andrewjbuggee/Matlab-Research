%% Simple viewing of Land-Ocean Masks that could be used for my retrieval

% Andrew John Buggee

%% Using the IMERG land sea mask (https://gpm.nasa.gov/data/directory/imerg-land-sea-mask-netcdf)

filename = 'IMERG_land_sea_mask.nc';

info = ncinfo(filename);

% read in lat long and mask
lat = ncread('IMERG_land_sea_mask.nc', "lat");
lon = ncread('IMERG_land_sea_mask.nc', "lon");
landSeaCover = ncread('IMERG_land_sea_mask.nc','landseamask');

% create geoshow plot
[lat_mat, lon_mat] = meshgrid(lat, lon);

% set up the type of global projection
axesm('apianus', 'Frame', 'on', 'Grid','on');
geoshow(lat_mat, lon_mat, landSeaCover);
cb = colorbar('southoutside');
set(gcf,'Position', [0 0 1400 1400])

% try imagesc
figure; imagesc(lat, lon, landSeaCover); colorbar
xlabel('Latitude'); ylabel('Longitude')
set(gcf,'Position', [0 0 1400 1400])



%% Read and View MODIS Land mask prodcut (https://lpdaac.usgs.gov/products/mcd12q1v006/


filename = 'MCD12Q1.A2020001.h05v10.006.2021359022853.hdf';

info = hdfinfo(filename);
