%% load a subset of suitable pixels and save the indices used



% By Andrew J. Buggee
%%

function pixels2use = subset_suitablePixels(inputs,modis)

suitablePixels_fileName = [inputs.savedCalculations_folderName,'suitablePixels.mat'];

folderName2Save = inputs.savedCalculations_folderName; % where to save the indices
numPixels2Calculate = inputs.pixels.num_2calculate; % number of pixels to use in our calcualtions

load(suitablePixels_fileName,'pixels');

if inputs.RT.use_custom_mie_calcs==false
    % ---------- Only use Pixels within mie_interpolate range ------
    % These ranges are listed in the pre-computed mie table from
    % libradtran.org
    re_mask = modis.cloud.effRadius17(pixels.res1km.index)>1 & modis.cloud.effRadius17(pixels.res1km.index)<25;
    tau_mask  = modis.cloud.optThickness17(pixels.res1km.index)>1 & modis.cloud.optThickness17(pixels.res1km.index)<80;
else
    error([newline, 'I dont know what range to use for re',newline])
end

% ---------- Only use Pixels with <10% uncertainty in retrieval ------
re_uncert_mask = modis.cloud.effRad_uncert_17(pixels.res1km.index)<10;
tau_uncert_mask = modis.cloud.optThickness_uncert_17(pixels.res1km.index)<10;


% ---------- Apply any additional tau_c constraints ------
tau_limit_mask = modis.cloud.optThickness17(pixels.res1km.index)>=inputs.pixels.tau_min_threshold & modis.cloud.optThickness17(pixels.res1km.index)<=inputs.pixels.tau_max_threshold;


% ---------- Apply any additional re constraints ------
re_limit_mask = modis.cloud.effRadius17(pixels.res1km.index)>=inputs.pixels.re_min_threshold & modis.cloud.effRadius17(pixels.res1km.index)<=inputs.pixels.re_max_threshold;



combined_mask = re_mask & tau_mask & re_uncert_mask & tau_uncert_mask & tau_limit_mask & re_limit_mask;

% % find the index value of all the ones!
% index = find(combined_mask);
% 
% % convert this to the indicies we started with. Which pixels that were
% % found in suitable_pixels also meet to requirements above?
% index_suit_pix = pixels.res1km.index(index);

index_suit_pix = pixels.res1km.index(combined_mask);


[row, col] = ind2sub(pixels.res1km.size, index_suit_pix);

% ------------------------------------------------------------
pixels_available.res1km.index = index_suit_pix;
pixels_available.res1km.row = row;
pixels_available.res1km.col = col;


% for all calculations, we need the 1km resolution pixel locations. So that
% is what we will use. We only use the 500 meter resolution pixel locaiton
% to show which pixels in the EV data set we used.

numSuitablePixels = length(pixels_available.res1km.index);


if numSuitablePixels > numPixels2Calculate
    
    rand_indices = randi([1, numSuitablePixels],numPixels2Calculate,1); % generate random numbers to choose pixels
    
    pixels2use.res1km.row = pixels_available.res1km.row(rand_indices); % row positions
    pixels2use.res1km.col = pixels_available.res1km.col(rand_indices); % column positions
    
    

    pixels2use.res1km.size = [size(modis.EV1km.radiance(:,:,1),1), size(modis.EV1km.radiance(:,:,1),2)];
    %pixels2use.res500m.size = pixels.res500m.size;
    
    % lets map 1 km pixels to 500 meter pixels location

    %[pixels2use.res500m.row, pixels2use.res500m.col] = cartesian_mapping_1000_to_500_pixels(pixels2use);


    
elseif numSuitablePixels < numPixels2Calculate
    
    % if we ask for more than there are, just use them all
    
    disp(['Number of Suitable Pixels is less than what you asked for.',...
        'Writing INP files for all suitabl epixels...']);
    
    pixels2use = pixels_available;
    

    
    
elseif numSuitablePixels == numPixels2Calculate
    
    pixels2use = pixels_available;
    
    
    
    
else
    
    error('Something is wrong with the number of pixels desired')
    
end


% --- Load the geometry settings for each pixel ---

pixel_rows = pixels2use.res1km.row;
pixel_cols = pixels2use.res1km.col;

% save the index
pixels2use.res1km.index = sub2ind([size(modis.EV1km.reflectance,1), size(modis.EV1km.reflectance,2)], pixel_rows, pixel_cols);

for ii = 1:length(pixel_rows)

    pixels2use.res1km.geometry.sza(ii) = modis.solar.zenith(pixel_rows(ii),pixel_cols(ii));
    pixels2use.res1km.geometry.saz(ii) = modis.solar.azimuth(pixel_rows(ii),pixel_cols(ii));
    
    % we need the cosine of the zenith viewing angle
    pixels2use.res1km.geometry.umu(ii) = round(cosd(double(modis.sensor.zenith(pixel_rows(ii),pixel_cols(ii)))),3); % values are in degrees
    pixels2use.res1km.geometry.phi(ii) = modis.sensor.azimuth(pixel_rows(ii),pixel_cols(ii));
    
end


% Save pixels2use and inputs
save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')




end