function f = plot_R17_and_T17_withUnertainty(modis,inputs)


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

% Plot locations where these thresholds are true

radiance_band1 = modis.EV1km.radiance(:,:,1);
reflectance_band1 = modis.EV1km.reflectance(:,:,1);

% lets look at the radiance for liquid water clouds only
radiance_band1_liquidWater = radiance_band1(liquidWater_mask);
lat_liquidWater = modis.geo.lat(liquidWater_mask);
long_liquidWater = modis.geo.long(liquidWater_mask);


f = figure; 
geoscatter(lat_liquidWater, long_liquidWater, 5, radiance_band1_liquidWater, '.')
title('Liquid Water Cloud Detected - Radiance 650 \mum band')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Radiance (W/m^{2}/\mum/sr)')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1500 700])



% lets look at the reflectance for the liquid water pixels with a certain
% optical depth threshold
reflectance_band1_combinedMask = reflectance_band1(combined_mask);
lat_combinedMask = modis.geo.lat(combined_mask);
long_combinedMask = modis.geo.long(combined_mask);


f = figure;subplot(1,2,1)
geoscatter(lat_combinedMask, long_combinedMask, 5, reflectance_band1_combinedMask, '.')
title(['Liquid Water Clouds w/ \tau_c>',num2str(tauThreshold),' - Reflectance 650 \mum band'])
cb = colorbar;
set(get(cb, 'label'), 'string', 'Reflectance (sr^{-1})')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')


% lets look at the reflectance for the liquid water pixels with a certain
% optical depth threshold
reflectance_band7 = modis.EV1km.reflectance(:,:,7);
reflectance_band7_combinedMask = reflectance_band7(combined_mask);


subplot(1,2,2)
geoscatter(lat_combinedMask, long_combinedMask, 5, reflectance_band7_combinedMask, '.')
title(['Liquid Water Clouds w/ \tau_c>',num2str(tauThreshold),' - Reflectance 2130 \mum band'])
cb = colorbar;
set(get(cb, 'label'), 'string', 'Reflectance (sr^{-1})')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1800 700])




% Lets look at the optical thickness estimates for all liquid water pixels
f = figure; subplot(1,2,1)
geoscatter(lat_combinedMask, long_combinedMask, 5, modis.cloud.optThickness17(combined_mask), '.')
title('MODIS Optical Thickness Retrieval')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Optical Thickness')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')



% Lets look at the droplet radius retrieval estimates for all liquid water pixels
subplot(1,2,2)
geoscatter(lat_combinedMask, long_combinedMask, 5, modis.cloud.effRadius17(combined_mask), '.')
title('MODIS Effective Radius Retrieval')
cb = colorbar;
set(get(cb, 'label'), 'string', 'Effective Radius (\mum)')
set(gca, 'FontSize',18)
set(gca, 'FontWeight', 'bold')
set(f, 'Position', [0 0 1800 700])


% --------------------------------------------
% ---------- DEFINE THE ERRORBARS ------------
% --------------------------------------------

re17 = modis.cloud.effRadius17(combined_mask);
T17 = modis.cloud.optThickness17(combined_mask);

xneg = re17.*(modis.cloud.effRad_uncert_17(combined_mask)./100);
xpos = xneg;
yneg = T17.*(modis.cloud.optThickness_uncert_17(combined_mask)./100);
ypos = xneg;

% ----- THERE ARE TWO MANY PIXELS!! FIND A RANDOM SAMPLE -----

colormap parula
num2plot = 1e3;
% 
% if numel(re17)>num2plot
%     
%     
%     rand_ind = randi([1,numel(re17)],num2plot,1);
%     
%     f = figure; errorbar(re17(rand_ind),T17(rand_ind),yneg(rand_ind),ypos(rand_ind),xneg(rand_ind),xpos(rand_ind),'o')
%     xlabel('r_{e}^{16} (\mu m)')
%     ylabel('\tau^{17}')
%     title('MODIS r_{e} and \tau_{c} Estimates')
%     grid on; grid minor
%     set(f, 'Position', [0 0 1000 400])
%     set(gca, 'YDir', 'reverse')
%     
% else
%     f = figure; errorbar(re17,T17,yneg,ypos,xneg,xpos,'o')
%     xlabel('r_{e}^{17} (\mu m)')
%     ylabel('\tau^{17}')
%     title('MODIS r_{e} and \tau_{c} Estimates')
%     grid on; grid minor
%     set(f, 'Position', [0 0 1000 400])
%     set(gca, 'YDir', 'reverse')
%     
% end



end

