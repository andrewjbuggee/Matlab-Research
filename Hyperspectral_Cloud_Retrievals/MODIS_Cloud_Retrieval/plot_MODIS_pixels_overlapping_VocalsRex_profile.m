%% Plot the VOCALS-REX flight path and the overlapping MODIS pixels

% This will plot all of the MODIS pixels found using modisIndex_min_dist
% These are the unique pixels found closest to the VOCALS-REx profile.


% INPUTS:

% (1) MODIS data structure
% (2) VOCALS-REx data structure
% (3) MODIS folder name, to provide the date the data was recorded
% (4) retrieved_variable - name of the retrieved variable to use as the
% colormap of the MODIS pixels

% 'radius' - plots each MODIS pixel with a color associatied with the
% magnitude of effective radius retrieved at that pixel by MODIS

% 'optical_depth' - plots each MODIS pixel with a color associatied with the
% magnitude of cloud optical depth retrieved at that pixel by MODIS


% By Andrew John Buggee

%%

function plot_MODIS_pixels_overlapping_VocalsRex_profile(modis,vocalsRex, modisFolder, modisInputs, retrieved_variable)



figure;


% load the across and along pixel growth curve fits
if strcmp(whatComputer, 'anbu8374')==true

    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/',...
        'along_across_pixel_growth.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/along_across_pixel_growth.mat'])

end

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


% ---- What retrieved variable should we plot? -----
if strcmp(retrieved_variable, 'radius')==true
    % grab the indexes of the ordered set of effective radii for the data
    % being used
    data2plot = modis.cloud.effRadius17(vocalsRex.modisIndex_minDist);

elseif strcmp(retrieved_variable, 'optical_depth')==true
    % grab the indexes of the ordered set of optical depths for the data
    % being used
    data2plot = modis.cloud.optThickness17(vocalsRex.modisIndex_minDist);

else

    error([newline, 'I dont recognize the retrieved variable!', newline])
end

% Set NaNs to 0
data2plot(isnan(data2plot))=0;

[~, index_sort] = sort(reshape(data2plot, [],1), 'ascend');

% Check to see if all the data is unique
[data2plot_unique, ~, idx_unique] = unique(data2plot(index_sort));

% We span the colormap only by the number of unique data values
C = parula(length(data2plot_unique));

% And now let's repeat the same colors for the redundant data points
C = C(idx_unique,:);


for ii = 1:length(data2plot)

    % the distance to the corners depends on the viewing zenith angle
    along_length = along_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));
    across_length = across_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));

    % compute the distance from the center to each corner
    dist_to_corner = sqrt((along_length/2)^2 + (across_length/2)^2) * 1e3;  % meters


    % The azimuth angle in the reckon function is with respect to the local
    % meridian (north). So we have to add the inclination angle to the azimuth
    % angles listed above.

    [lat_corner, long_corner] = reckon(modis.geo.lat(vocalsRex.modisIndex_minDist(ii)),...
        modis.geo.long(vocalsRex.modisIndex_minDist(ii)),linspace(dist_to_corner, dist_to_corner,5),...
        azimuth_2_corners+terra_inclination, wgs84Ellipsoid);

    % Create a geopolyshape
    modis_polyshape = geopolyshape(lat_corner, long_corner);

    gp = geoplot(modis_polyshape);
    gp.EdgeAlpha = 0;
    % The color corresponds to the linear index_sort
    gp.FaceColor = C(ii==index_sort,:);
    gp.FaceAlpha = 0.7;



    hold on



end


% ----- Plot the Vocals-Rex data ------
if modisInputs.flags.useAdvection==true

    geoscatter(vocalsRex.lat_withAdvection, vocalsRex.long_withAdvection, 10, "red",'*')

else

    geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10, "red",'*')

end


cb = colorbar;

% ---- What retrieved variable should we plot? -----
if strcmp(retrieved_variable, 'radius')==true
    % grab the indexes of the ordered set of effective radii for the data
    % being used
    clim([min(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all'),...
        max(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all')])

    set(get(cb, 'label'), 'string', '$r_e \; (\mu m)$','Interpreter','latex', 'Fontsize',28)


    % load the across and along pixel growth curve fits
    if strcmp(whatComputer, 'anbu8374')==true

        title(['MODIS Effective Radius - ', modisFolder(97:end-1)],'Interpreter','latex', 'FontSize', 30)

    elseif strcmp(whatComputer, 'andrewbuggee')==true

        title(['MODIS Effective Radius - ', modisFolder(113:end-1)],'Interpreter','latex', 'FontSize', 30)

    end



elseif strcmp(retrieved_variable, 'optical_depth')==true
    % grab the indexes of the ordered set of optical depths for the data
    % being used
    clim([min(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist), [], 'all'),...
        max(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist), [], 'all')])

    set(get(cb, 'label'), 'string', '$\tau_c$','Interpreter','latex', 'Fontsize',28)


    % load the across and along pixel growth curve fits
    if strcmp(whatComputer, 'anbu8374')==true

        title(['MODIS Optical Depth - ', modisFolder(97:end-1)],'Interpreter','latex', 'FontSize', 30)

    elseif strcmp(whatComputer, 'andrewbuggee')==true

        title(['MODIS Optical Depth - ', modisFolder(113:end-1)],'Interpreter','latex', 'FontSize', 30)

    end


else

    error([newline, 'I dont recognize the retrieved variable!', newline])
end



set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 600 600])


if modisInputs.flags.useAdvection==true
% Create textbox
annotation('textbox',[0.0593333333333334 0.0266666666666667 0.220666666666667 0.0583333333333333],...
    'String','With Advection',...
    'LineWidth',2,...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

elseif modisInputs.flags.useAdvection==false
% Create textbox
annotation('textbox',[0.0593333333333334 0.0266666666666667 0.267333333333333 0.0583333333333333],...
    'String','Without Advection',...
    'LineWidth',2,...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

end




end