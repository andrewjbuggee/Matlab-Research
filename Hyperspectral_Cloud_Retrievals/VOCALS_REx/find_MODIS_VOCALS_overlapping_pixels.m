function pixels2use = find_MODIS_VOCALS_overlapping_pixels(modis, inputs, vocalsRex)

% Define folder to save calculations
folderName2Save = inputs.savedCalculations_folderName; % where to save the indices

% NOT USING SUITABLE PIXELS FOR NOW

% ----------------------------------------------------------------
% ------- FIND PIXELS THAT MEET REQUIREMENTS LISTED BELOW --------
% ----------------------------------------------------------------
 % 2 is the value designated for liquid water
liquidWater_mask = modis.cloud.phase == 2;

% create tau mask based on threshold
tauThreshold = inputs.pixels.tau_min_threshold;
%tauThreshold = 0;

% finds clouds with an optical thickness of a certain value and an
% uncertainty less than the definition below
uncertaintyLimit = 10;                              % percentage


% Let's check to see that the MODIS data that coincides with the cloud
% profile chosen within vocals-Rex matches these requirements
indexes2check = vocalsRex.modisIndex_minDist;

for nn = 1:numel(indexes2check)


    if liquidWater_mask(indexes2check(nn))==true && modis.cloud.optThickness17(indexes2check(nn))>=tauThreshold && modis.cloud.effRad_uncert_17(indexes2check(nn))<uncertaintyLimit

    else
        error([newline,'The MODIS pixel chosen for comparison with the VOCALS-REx measurement is incompatible',newline,...
            'Either this pixel didnt detect liquid water, or it didnt surpass the tau threshold',newline])

    end

end


% ------------------------------------------------------------------------
% ----------------------------- OLD CODE ---------------------------------
% ------------------------------------------------------------------------


% % Now we step through each vocals Rex location and find the pixel in the
% % MODIS scene that is closest to it
% 
% num_vocals_locations = length(vocalsRex.latitude(index));
% 
% 
% % only use the indexs where the plane is rising into a cloud
% vocalsLat = vocalsRex.latitude(index);
% vocalsLong = vocalsRex.longitude(index);
% 
% index_closest_pixel = zeros(1,num_vocals_locations);
% 
% for nn = 1:num_vocals_locations
% 
%     dist = sqrt((modis.geo.lat - vocalsLat(nn)).^2 + (modis.geo.long - vocalsLong(nn)).^2);
%     [~,index_closest_pixel(nn)] = min(dist,[],"all");
%     
% end
% 
% % There are bound to be redundant values. Take only the unique indexes
% indexes2keep = unique(index_closest_pixel);
% 
% % Now we need to sort through and make sure to keep only the indexes with a
% % liquid water phase, a tau above our threshold and the retrieval below our
% % uncertainty limit
% index_threshold = modis.cloud.phase(indexes2keep)==2 & modis.cloud.optThickness17(indexes2keep)>=tauThreshold & modis.cloud.effRad_uncert_17(indexes2keep)<=uncertaintyLimit;
% 
% % This is the final list! Keep all of these pixels!


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------



% store the linear index, the rows and the columns
pixels2use.res1km.linearIndex = indexes2check;

[pixels2use.res1km.row, pixels2use.res1km.col] = ind2sub(size(modis.geo.lat), pixels2use.res1km.linearIndex);
    



% --- Load the geometry settings for each pixel ---

pixels2use.res1km.geometry.sza = modis.solar.zenith(pixels2use.res1km.linearIndex);
pixels2use.res1km.geometry.saz = modis.solar.azimuth(pixels2use.res1km.linearIndex);
pixels2use.res1km.geometry.umu = round(cosd(double(modis.sensor.zenith(pixels2use.res1km.linearIndex))),3); % values are in degrees
pixels2use.res1km.geometry.phi = modis.sensor.azimuth(pixels2use.res1km.linearIndex);

% Save the pixels to a file, and save the geometry in the pixels2use
% structure

save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')








end


