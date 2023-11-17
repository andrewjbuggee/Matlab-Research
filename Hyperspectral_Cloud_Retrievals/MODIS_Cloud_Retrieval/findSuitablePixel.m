%% ----- Algorithm to search for a suitable Pixel -----

% This function looks through a modis data set to find a pixel that has a
% thermodynamic phase of liquid water, a certain optical depth threshold,
% and is surrounded by many other pixels like it - thus it is not an edge case
% For now, the pixel MUST be over ocean.

% save both 1km resolution pixels and we interpolate to find how these map
% to 500 meter resolution pixels


% By Andrew J. Buggee

%%

function [pixels] = findSuitablePixel(modis,inputs)

% ---------------------------------------------------------------------
% LETS DEFINE MINIMUM VALUES FOR THE OPTICAL DEPTH AND THE UNCERTAINTY
% PIXELS ABOVE THESE MINIMUM VALUES WILL BE USED IN THE TBLUT ALGORITHM
% ---------------------------------------------------------------------


% find pixels above a certain optical depth
tau_min_threshold = inputs.pixels.tau_min_threshold;

% find pixels below a certain optical depth
tau_max_threshold = inputs.pixels.tau_max_threshold;

% find pixels above a certain effective radius
re_min_threshold = inputs.pixels.re_min_threshold;

% find pixels below a certain effective radius
re_max_threshold = inputs.pixels.re_max_threshold;


% create logical mask for phase

% find pixels that detect liquid water
liquidWater_mask = modis.cloud.phase == 2; % 2 is the value designated for liquid water


% create tau mask based on threshold

% finds clouds with an optical thickness of a certain value and an
% uncertainty less than the definition below.
%
% Uncertainties are in percentages
uncertaintyLimit = 10;                              % percentage

tau_mask = modis.cloud.optThickness17 >= tau_min_threshold & modis.cloud.optThickness17 <= tau_max_threshold & modis.cloud.optThickness_uncert_17<uncertaintyLimit;


% Find pixels with an effective radius of at least 0 and an uncertainty
% less than the amount defined below. Uncertainties are in percentages
uncertaintyLimit = 10;
re_mask = modis.cloud.effRadius17>=re_min_threshold & modis.cloud.effRadius17<=re_max_threshold & modis.cloud.effRad_uncert_17<uncertaintyLimit;        % find values greater than 0

% find where there is overlap

combined_mask = logical(liquidWater_mask .* tau_mask.*re_mask);
% --------------------------------------------------------------


% ---- Define distance between liquid and non-liquid pixels ----
% set up a distance threshold, which is the number of pixels from a boundry
% where the phase changes

dist_threshold = 2; % this defines the closest possible distance to a 0 pixel in the combined mask


% ----- Remove edge pixels (WHY?!) -----
% we don't want any edge pixels. So the pixels also have to be atleast border_threshold away from the border of the image
% This will reduce the number of computations! To do this we will turn the
% 1's around the border into 0's

border_threshold = dist_threshold; % distance away from image edge a pixel has to be in order to be considered

border_mask = ones(size(combined_mask));
border_mask([1:border_threshold,(size(border_mask,1)-border_threshold+1):end],:) = 0; % convert top and bottom border into zeros
border_mask(:,[1:border_threshold,(size(border_mask,2)-border_threshold+1):end],:) = 0; % convert left and right border into zeros

combined_mask = combined_mask .* border_mask;
% --------------------------------------------------------------


% ---- Only select pixels over ocean -----

% test locations from combined_mask
coastal_res = 20;    % 1 is low resolution, 10 is decently high resolution
make_plot = 0;  %0 = no plot, 1 = plot
isOcean = land_or_ocean(double(modis.geo.lat(:)),double(modis.geo.long(:)),...
    coastal_res,make_plot);
% reshape ocean so it's the same size as combined_mask
isOcean = reshape(isOcean, size(modis.geo.lat,1), size(modis.geo.long,2));

% create new combined mask so that we only use pixels over ocean
combined_mask = combined_mask .* isOcean;
% --------------------------------------------------------------


% ---- Check 8 neighboring pixels for different phase -----

% this will check to see if the 8 neighboring pixels are all liquid water.
% If they are, we keep the pixel. Otherwise, we don't keep it.


% Find all the pixels the meet our requirements

index_ones = find(combined_mask); % finds the indices of the non-zero elements - these are the pixels that meet all the above requirements
%index_zeros = find(~logical(combined_mask)); % finds the indices of the zero elements - these are the pixels that DONT meet the above requirements

%[row0,col0] = ind2sub(size(combined_mask),index_zeros); % convert to row column indices
[row1,col1] = ind2sub(size(combined_mask),index_ones); % conver to row column indices







% ----- Working version of the suitable pixels code -----

%indexes_2keep = [];
%
% while isempty(indexes_2keep) == true
% 
% 
%     % ---- PARALLEL FOR LOOP -----
%     for ii = 1:length(row1)
% 
% 
% 
%         dist2pixel = sqrt((row1(ii) - row0).^2 + (col1(ii) - col0).^2);
% 
%         if min(dist2pixel,[],'all') >= dist_threshold
% 
%             indexes_2keep = [indexes_2keep; index_ones(ii)];
%         end
% 
%         % ----------------------------------------------
% 
% 
%     end
% 
%     if isempty(indexes_2keep) == true
%         disp(['No Suitable pixels found for a distance threshold of ',num2str(dist_threshold)])
%         dist_threshold = floor(dist_threshold/2);
% 
%     elseif isempty(indexes_2keep) == false
%         disp([num2str(numel(indexes_2keep)),' suitable pixels found with a distance threshold of ',num2str(dist_threshold)])
%     end
% 
% 
% end





% ---------------------------------------------------------------
%     ----- Working on Faster code? ------
%     Let's try new code to speed up this process.
%     We need to check that all the requirements are met for all pixels
%     within the distance threshold.

% create an empty array for the indexes we will keep

indexes_2keep = [];

% The size of the array we will test will be always be the same
pixels2check_size = (2*dist_threshold + 1)^2;

% Start by stepping through the pixels that meet the above requirements
for ii = 1:length(row1)

    % set up a mesh grid cenetered on the pixel found to match all the
    % requirements above.

    % check the see if there are any 0's in the combined mask within the
    % threshold arond the pixel found to meet the above requirements
    pixels2check = combined_mask(row1(ii)-dist_threshold:row1(ii)+dist_threshold, col1(ii)-dist_threshold:col1(ii)+dist_threshold);

    % If all the pixels within our distance threshold meet the above constraints, we keep the index 
    if length(find(pixels2check))== pixels2check_size

        % keep this index!
        indexes_2keep = [indexes_2keep, index_ones(ii)];

    end
    
end





% check to see if we found any indexes to keep
if isempty(indexes_2keep)==true

    error([newline, 'No suitable pixels found!', newline])

else
    disp([newline, num2str(numel(indexes_2keep)),' suitable pixels found with a distance threshold of ',num2str(dist_threshold), newline])


end


%     -----------------------------------------------------------
%     -----------------------------------------------------------


% lets keep both the 1km resolution pixel location, and interpolate to get pixel
% locations for 500 meter resolution data

[pixels.res1km.row, pixels.res1km.col] = ind2sub(size(combined_mask),indexes_2keep);

% save the indices
pixels.res1km.index = indexes_2keep;
% save the size of each resolution swath
%pixels.res500m.size = size(modis.EV500m.reflectance(:,:,1)); % the 500 meter resolution image swath
pixels.res1km.size = size(modis.EV1km.reflectance(:,:,1)); % the 1km resolution image swath

% lets map 1 km pixels to 500 meter pixels location

%[pixels.res500m.row, pixels.res500m.col] = cartesian_mapping_1000_to_500_pixels(pixels);

end






