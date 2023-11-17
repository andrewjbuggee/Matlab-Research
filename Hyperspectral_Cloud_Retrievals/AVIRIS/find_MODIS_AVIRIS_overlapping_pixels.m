function pixels2use = find_MODIS_AVIRIS_overlapping_pixels(modis, inputs, aviris)

% Define folder to save calculations
folderName2Save = inputs.savedCalculations_folderName; % where to save the indices

% NOT USING SUITABLE PIXELS FOR NOW

% ----------------------------------------------------------------
% ------- FIND PIXELS THAT MEET REQUIREMENTS LISTED BELOW --------
% ----------------------------------------------------------------
 % 2 is the value designated for liquid water
liquidWater_mask = modis.cloud.phase == 2;

% create tau mask based on threshold
tauThreshold = inputs.pixels.tauThreshold;
%tauThreshold = 0;

% finds clouds with an optical thickness of a certain value and an
% uncertainty less than the definition below
uncertaintyLimit = 10;                              % percentage

tau_mask = modis.cloud.optThickness17 >= tauThreshold & modis.cloud.optThickness_uncert_17<uncertaintyLimit;



% Find pixels with am effective radius of at least 0 and an uncertainty
% less than the amount defined below
uncertaintyLimit = 10;
re_mask = modis.cloud.effRadius17>=0 & modis.cloud.optThickness_uncert_17<uncertaintyLimit;        % find values greater than 0

% find where there is overlap

combined_mask = logical(liquidWater_mask .* tau_mask.*re_mask);

suitablePixel_lat = modis.geo.lat(combined_mask);
suitablePixel_long = modis.geo.long(combined_mask);

% ----------------------------------------------------------------
% ----------------------------------------------------------------

% We need the linear index of all the suitable pixels found above. This is
% how we map back to the pixels of interest
linearIndex = find(combined_mask);


% I want to find MODIS pixels that are within the AVIRIS scene.
% Lets start be defining the AVIRIS scene border


% AVIRIS defines longitude between -180 and 180 degrees, where values less
% than 0 are west of the Greenwich merdian. MODIs defines lon
avirisLatBorder = [aviris.position(:,1,2),aviris.position(:,end,2)];
avirisLongBorder = [aviris.position(:,1,1),aviris.position(:,end,1)];


% Lets step through each latitude and find where the data is bounded by
% AVIRIS data in longitude

dLat_MODIS = 0.01;           % This should be broad enough to find

% break up the aviris border into N sections eqaul to the total length
% divided by the MODIS step size
N = floor((avirisLatBorder(1,1) - avirisLatBorder(end,1))/dLat_MODIS);

% Lets create the latitude vector that we will step through
lat_vector = avirisLatBorder(end,1):dLat_MODIS:avirisLatBorder(1,1);

% create the longitude vector that corresponds to the latitude vector
% do this by finding the longitudinal midpoint between the latitude
% boundaires. Because we're searching within latitude bands, the longitude
% vector is 1 less in length than the latitude vector. The longtidue vector
% represents the layers inbetween lines of constant latitude.

long_vectorLeft = linspace(avirisLongBorder(end,1), avirisLongBorder(1,1),N);

% shift the longitude vector to obtain the midpoint of each layer
long_vectorLeft = long_vectorLeft + (long_vectorLeft(2)-long_vectorLeft(1))/2;

% Do the same thing for the right hand side
long_vectorRight = linspace(avirisLongBorder(end,2), avirisLongBorder(1,2),N);

% shift the longitude vector to obtain the midpoint of each layer
long_vectorRight = long_vectorRight + (long_vectorRight(2)-long_vectorRight(1))/2;


% start an empty vector which will contain the locations of the selected
% pixels
linearIndex_selectedPixels = [];

for ii = 1:(N-1)
    
    % Step through each latitude, and find MODIS pixels that lie between
    % the MODIS longitude barriers
    
    index_lat = suitablePixel_lat>=lat_vector(ii) & suitablePixel_lat<lat_vector(ii+1);
    
    % Now find the pixels that are within the longitude boundaries
    
    index_long = suitablePixel_long>=long_vectorLeft(ii) & suitablePixel_long<=long_vectorRight(ii);
    
    combined_index = logical(index_lat.*index_long);
    
    % Now grab the linear indices of the pixels within the AVIRIS scene
    linearIndex_selectedPixels = [linearIndex_selectedPixels;  linearIndex(combined_index)];
    
    
end

% if we didn't find any pixels, declare an error so the user knows
if isempty(linearIndex_selectedPixels)==true
    error([newline, 'No MODIS pixels with the minimum specifications defined were found within ',...
        'the AVIRIS image swath!', newline]);
end

% We only want to calculate the number of pixels that was specified in the
% inputs

if inputs.pixels.num_2calculate>=length(linearIndex_selectedPixels)
    
    % store the linear index, the rows and the columns
    pixels2use.res1km.linearIndex = linearIndex_selectedPixels;
    [pixels2use.res1km.row, pixels2use.res1km.col] = ind2sub(size(modis.geo.lat), linearIndex_selectedPixels);
    
else
    
    % if there are more pixels found than is desired for calculation, take
    % a random subset
    random_set = randi([1, length(linearIndex_selectedPixels)], inputs.pixels.num_2calculate,1);
    
    % store the linear index, the rows and the columns
    pixels2use.res1km.linearIndex = linearIndex_selectedPixels(random_set);
    [pixels2use.res1km.row, pixels2use.res1km.col] = ind2sub(size(modis.geo.lat), linearIndex_selectedPixels(random_set));
    
    
    
    
    
    
end



% --- Load the geometry settings for each pixel ---

pixel_rows = pixels2use.res1km.row;
pixel_cols = pixels2use.res1km.col;

for ii = 1:length(pixel_rows)
    
    pixels2use.res1km.geometry.sza(ii) = modis.solar.zenith(pixel_rows(ii),pixel_cols(ii));
    pixels2use.res1km.geometry.saz(ii) = modis.solar.azimuth(pixel_rows(ii),pixel_cols(ii));
    
    % we need the cosine of the zenith viewing angle
    pixels2use.res1km.geometry.umu(ii) = round(cosd(double(modis.sensor.zenith(pixel_rows(ii),pixel_cols(ii)))),3); % values are in degrees
    pixels2use.res1km.geometry.phi(ii) = modis.sensor.azimuth(pixel_rows(ii),pixel_cols(ii));
    
end

% Save the pixels to a file, and save the geometry in the pixels2use
% structure

save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')




end


