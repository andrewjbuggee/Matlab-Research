%% Convert MODIS Lat/Lon coordinates to a geo-polygon shape


% By Andrew John Buggee

%%

lat = modis.geo.lat(vocalsRex.modisIndex_minDist(1));
long = modis.geo.long(vocalsRex.modisIndex_minDist(1));

% Assume these are the center coordinates of the pixel
% First, create a perfect square assuming, when MODIS is looking nadir,
% that the pixel is 1km by 1 km on the ground.
% Thus, the diagonal from the center to any corner is 1/sqrt(2) km long
dist_to_corners = 1/sqrt(2);        % km

% Using the 'reckon' function, we con move to each corner from the center
% position. But we need to know the orientation of the MODIS instrument
% with respect to the meriodional lines.

terra_inclination = 98.2098;        % degrees
aqua_inclination = 	98.1987;        % degrees

azimuth_2_corners = [45, 135, 225, 315, 45];

% The azimuth angle in the reckon function is with respect to the local
% meridian (north). So we have to add the inclination angle to the azimuth
% angles listed above.

[lat_corner, long_corner] = reckon(lat, long, linspace(dist_to_corners, dist_to_corners,5)*1e3,...
    azimuth_2_corners+terra_inclination, wgs84Ellipsoid);

% Create a geopolyshape
modis_polyshape = geopolyshape(lat_corner, long_corner);

%% Create geoplot

figure;
geoplot(modis_polyshape)



%% Plot many modis pixels at once

% load the across and along pixel growth curve fits
load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
    'MODIS_Cloud_Retrieval/along_across_pixel_growth.mat'])

% Assume these are the center coordinates of the pixel
% First, create a perfect square assuming, when MODIS is looking nadir,
% that the pixel is 1km by 1 km on the ground.
% Thus, the diagonal from the center to any corner is 1/sqrt(2) km long
%dist_to_corners = 1/sqrt(2);        % km

% Using the 'reckon' function, we con move to each corner from the center
% position. But we need to know the orientation of the MODIS instrument
% with respect to the meriodional lines.

terra_inclination = 98.2098;        % degrees
aqua_inclination = 	98.1987;        % degrees

azimuth_2_corners = [45, 135, 225, 315, 45];

% Create the Color pallete
% n_rows = size(modis.geo.lat,1)/10;
%n_cols = size(modis.geo.lat,2);
row_start = 1500;
row_end = 2000;
col_start = 1;
col_end = 100;


% grab the indexes of the ordered set of effective radius for the data
% being used
data2plot = modis.cloud.effRadius17(row_start:row_end, col_start:col_end);
% Set NaNs to 0
data2plot(isnan(data2plot))=0;

[~, index_sort] = sort(reshape(data2plot, [],1), 'ascend');
%[r,c] = ind2sub([n_rows, n_cols], index_sort);
C = parula(length(index_sort));

figure;

for rr = row_start:row_end
    for cc = col_start:col_end

        % the distance to the corners depends on the viewing zenith angle
        along_length = along_scan(modis.sensor.zenith(rr,cc));
        across_length = across_scan(modis.sensor.zenith(rr,cc));

        % compute the distance from the center to each corner
        dist_to_corner = sqrt((along_length/2)^2 + (across_length/2)^2) * 1e3;  % meters


        % The azimuth angle in the reckon function is with respect to the local
        % meridian (north). So we have to add the inclination angle to the azimuth
        % angles listed above.

        [lat_corner, long_corner] = reckon(modis.geo.lat(rr,cc), modis.geo.long(rr,cc),...
            linspace(dist_to_corner, dist_to_corner,5), azimuth_2_corners+terra_inclination, wgs84Ellipsoid);

        % Create a geopolyshape
        modis_polyshape = geopolyshape(lat_corner, long_corner);

        gp = geoplot(modis_polyshape);
        gp.EdgeAlpha = 0;
        % The color corresponds to the linear index_sort
        idx_lin = sub2ind(size(data2plot), rr - row_start +1, cc - col_start + 1);
        gp.FaceColor = C(idx_lin==index_sort,:);
        gpFaceAlpha = 0.9;

        

        hold on


    end
end

cb = colorbar;
%clim([min(modis.cloud.effRadius17(1:n_rows, 1:n_cols)), max(modis.cloud.effRadius17(1:n_rows, 1:n_cols))])
clim([0, max(modis.cloud.effRadius17(row_start:row_end, col_start:col_end), [], 'all')])

set(get(cb, 'label'), 'string', '$r_e \; (\mu m)$','Interpreter','latex', 'Fontsize',22)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
