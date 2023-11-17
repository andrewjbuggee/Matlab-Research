function f = plot_MODIS_re_and_Tau_with_VOCALS_overlap(modis,inputs, vocalsRex)


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


% Plot locations where these thresholds are true

radiance_band1 = modis.EV1km.radiance(:,:,1);
reflectance_band1 = modis.EV1km.reflectance(:,:,1);

% lets look at the radiance for liquid water clouds only
radiance_band1_liquidWater = radiance_band1(liquidWater_mask);
lat_liquidWater = modis.geo.lat(liquidWater_mask);
long_liquidWater = modis.geo.long(liquidWater_mask);



% Extract just the border of the AVIRIS image
vocalsLatBorder = vocalsRex.latitude;
vocalsLongBorder = vocalsRex.longitude;

latLimits = [min(vocalsLatBorder)-2, max(vocalsLatBorder)+2];
longLimits = [min(vocalsLongBorder)-2, max(vocalsLongBorder)+2];

% Define the marker size
modisMarker = 40;
VocalsMarker = 80;


% Lets look at the optical thickness estimates for all liquid water pixels
f = figure; subplot(1,2,1)
geoscatter(lat_combinedMask, long_combinedMask, modisMarker, modis.cloud.optThickness17(combined_mask), '.')
title('MODIS Optical Thickness Retrieval','Interpreter','latex','FontSize',30)
cb = colorbar;
set(get(cb, 'label'), 'string', 'Optical Thickness','Interpreter','latex', 'Fontsize',24)
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
hold on;
geoscatter(vocalsLatBorder, vocalsLongBorder, VocalsMarker, 'r', '.')
geolimits(latLimits, longLimits)



% Lets look at the droplet radius retrieval estimates for all liquid water pixels
subplot(1,2,2)
geoscatter(lat_combinedMask, long_combinedMask, modisMarker, modis.cloud.effRadius17(combined_mask), '.')
title('MODIS Effective Radius Retrieval','Interpreter','latex', 'FontSize',30)
cb = colorbar;
set(get(cb, 'label'), 'string', 'Effective Radius ($\mu m$)','Interpreter','latex','Fontsize',24)
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1800 700])
hold on;
geoscatter(vocalsLatBorder, vocalsLongBorder, VocalsMarker, 'r', '.')
geolimits(latLimits, longLimits)
% Create textbox describing the AVIRIS data location
dim = [.45 .68 .3 .3];
str = 'VOCALS-REx Flight Path in red';
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','r',...
           'FontWeight','bold','FontSize',14, 'EdgeColor','k');
% Create textbox describing the optical depth threshold
dim = [.49 0.6 .3 .3];
str = ['\tau_{c} > ',num2str(tauThreshold)];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
           'FontWeight','bold','FontSize',14, 'EdgeColor','k');
       
       
       
% ---------------------------------------------         
% ---------------------------------------------       
% Lets also plot the cloud top height and phase
% ---------------------------------------------  
% ---------------------------------------------  

% % Lets look at the optical thickness estimates for all liquid water pixels
% f = figure; subplot(1,2,1)
% geoscatter(lat_combinedMask, long_combinedMask, modisMarker, modis.cloud.topHeight(combined_mask), '.')
% title('MODIS Cloud Top Height (m)')
% cb = colorbar;
% set(get(cb, 'label'), 'string', 'Cloud Top Height (m)')
% set(gca, 'FontSize',18)
% set(gca, 'FontWeight', 'bold')
% hold on;
% geoscatter(vocalsLatBorder, vocalsLongBorder, VocalsMarker, 'r', '.','MarkerFaceAlpha',0.1)
% geolimits(latLimits, longLimits)
% 
% 
% 
% % Lets look at the droplet radius retrieval estimates for all liquid water pixels
% subplot(1,2,2)
% geoscatter(modis.geo.lat(liquidWater_mask), modis.geo.long(liquidWater_mask), modisMarker, modis.cloud.phase(liquidWater_mask), '.')
% title('MODIS Detection of Liquid Water')
% 
% set(gca, 'FontSize',18)
% set(gca, 'FontWeight', 'bold')
% set(f, 'Position', [0 0 1800 700])
% hold on;
% geoscatter(vocalsLatBorder, vocalsLongBorder, VocalsMarker, 'r', '.','MarkerFaceAlpha',0.1)
% geolimits(latLimits, longLimits)
% % Create textbox describing the AVIRIS data location
% dim = [.45 .68 .3 .3];
% str = 'VOCALS-REx Flight Path in red';
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
%            'FontWeight','bold','FontSize',14, 'EdgeColor','k');
% % Create textbox describing the optical depth threshold
% dim = [.49 0.6 .3 .3];
% str = ['\tau_{c} > ',num2str(tauThreshold)];
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
%            'FontWeight','bold','FontSize',14, 'EdgeColor','k');
       
       




end
