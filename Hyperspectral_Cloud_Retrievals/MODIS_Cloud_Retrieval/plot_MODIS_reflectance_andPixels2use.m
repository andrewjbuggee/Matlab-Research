function f = plot_MODIS_reflectance_andPixels2use(modis,pixels2use)


liquidWater_mask = modis.cloud.phase == 2; % 2 is the value designated for liquid water


% Plot locations where these thresholds are true

%radiance_band1 = modis.EV1km.radiance(:,:,1);
reflectance_band1 = modis.EV1km.reflectance(:,:,1);

% lets look at the radiance for liquid water clouds only
%radiance_band1_liquidWater = radiance_band1(liquidWater_mask);
reflectance_band1_liquidWater = reflectance_band1(liquidWater_mask);

lat_liquidWater = modis.geo.lat(liquidWater_mask);
long_liquidWater = modis.geo.long(liquidWater_mask);


f = figure; 
geoscatter(lat_liquidWater, long_liquidWater, 5, reflectance_band1_liquidWater, '.')
title('Liquid Water Cloud Detected - Radiance 650 \mum band')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Reflectance (1/sr)')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1500 700])
hold on

% Plot the locations of the pixels used for calculation
for ii = 1:length(pixels2use.res1km.col)
    
    geoscatter(modis.geo.lat(pixels2use.res1km.row(ii),pixels2use.res1km.col(ii)), modis.geo.long(pixels2use.res1km.row(ii),pixels2use.res1km.col(ii)), 10, 'k*')
end
set(f, 'Position', [0 0 1500 700])

dim = [.6 .6 .3 .3];
str = '* show location of pixels used in TBLUT calculation';
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
           'FontWeight','bold','FontSize',14, 'EdgeColor','k');



end

