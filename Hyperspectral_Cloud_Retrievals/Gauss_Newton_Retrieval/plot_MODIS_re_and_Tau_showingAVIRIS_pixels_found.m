function f = plot_MODIS_re_and_Tau_showingAVIRIS_pixels_found(modis,inputs, aviris, pixels2use)


liquidWater_mask = modis.cloud.phase == 2; % 2 is the value designated for liquid water

% create tau mask based on threshold
tauThreshold = inputs.pixels.tauThreshold;
%tauThreshold = 5;

% finds clouds with an optical thickness of a certain value and an
% uncertainty less than the definition below
uncertaintyLimit = 10;                              % percentage

tau_mask = modis.cloud.optThickness17 >= tauThreshold & modis.cloud.optThickness_uncert_17<uncertaintyLimit;



% Find pixels with am effective radius of at least 0 and an uncertainty
% less than the amount defined below
uncertaintyLimit =10;
re_mask = modis.cloud.effRadius17>=0 & modis.cloud.optThickness_uncert_17<uncertaintyLimit;        % find values greater than 0

% find where there is overlap

combined_mask = logical(liquidWater_mask .* tau_mask.*re_mask);


% lets look at the reflectance for the liquid water pixels with a certain
% optical depth threshold
lat_combinedMask = modis.geo.lat(combined_mask);
long_combinedMask = modis.geo.long(combined_mask);



% Extract just the border of the AVIRIS image
avirisLatBorder = [aviris.position(:,1,2);aviris.position(:,end,2); aviris.position(1,:,2)'; aviris.position(end,:,2)'];
avirisLongBorder = [aviris.position(:,1,1);aviris.position(:,end,1); aviris.position(1,:,1)'; aviris.position(end,:,1)'];

latLimits = [min(avirisLatBorder)-0.1, max(avirisLatBorder)+0.1];
longLimits = [min(avirisLongBorder)-0.1, max(avirisLongBorder)+0.1];


% Lets look at the optical thickness estimates for all liquid water pixels
f = figure; subplot(1,2,1)
geoscatter(lat_combinedMask, long_combinedMask, 30, modis.cloud.optThickness17(combined_mask), '.')
title('MODIS Optical Thickness Retrieval')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Optical Thickness')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
hold on;
geoscatter(avirisLatBorder, avirisLongBorder, 1, 'r', '.','MarkerFaceAlpha',0.1)
geoscatter(modis.geo.lat(pixels2use.res1km.linearIndex), modis.geo.long(pixels2use.res1km.linearIndex), 30, 'k*')
geolimits(latLimits, longLimits)



% Lets look at the droplet radius retrieval estimates for all liquid water pixels
subplot(1,2,2)
geoscatter(lat_combinedMask, long_combinedMask, 30, modis.cloud.effRadius17(combined_mask), '.')
title('MODIS Effective Radius Retrieval')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Effective Radius (\mum)')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1800 700])
hold on;
geoscatter(avirisLatBorder, avirisLongBorder, 1, 'r', '.','MarkerFaceAlpha',0.1)
geoscatter(modis.geo.lat(pixels2use.res1km.linearIndex), modis.geo.long(pixels2use.res1km.linearIndex), 30, 'k*')
geolimits(latLimits, longLimits)
% Create textbox describing the AVIRIS data location
dim = [.15 .6 .3 .3];
str = 'AVIRIS data within red box';
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
           'FontWeight','bold','FontSize',14, 'EdgeColor','k');
% Create textbox describing the optical depth threshold
dim = [.49 0.6 .3 .3];
str = ['\tau_{c} > ',num2str(tauThreshold)];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
           'FontWeight','bold','FontSize',14, 'EdgeColor','k');



end
