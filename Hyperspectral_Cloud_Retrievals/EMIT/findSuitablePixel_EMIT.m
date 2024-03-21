%% ----- Algorithm to search for a suitable Pixels within EMIT data -----




% By Andrew J. Buggee

%%

function [pixels] = findSuitablePixel_EMIT(emit)






% ---- Only select pixels over ocean -----

% test locations from combined_mask
coastal_res = 20;    % 1 is low resolution, 10 is decently high resolution
make_plot = 0;  %0 = no plot, 1 = plot
isOcean = land_or_ocean(double(emit.lat(:)),double(emit.long(:)),...
    coastal_res,make_plot);
% reshape ocean so it's the same size as combined_mask
isOcean = reshape(isOcean, size(emit.lat,1), size(emit.long,2));

% create new combined mask so that we only use pixels over ocean
combined_mask = isOcean;
% --------------------------------------------------------------


% ---- Check 8 neighboring pixels for different phase -----

% this will check to see if the 8 neighboring pixels are all liquid water.
% If they are, we keep the pixel. Otherwise, we don't keep it.


% Find all the pixels the meet our requirements

pixels.index = find(combined_mask); % finds the indices of the non-zero elements - these are the pixels that meet all the above requirements
%index_zeros = find(~logical(combined_mask)); % finds the indices of the zero elements - these are the pixels that DONT meet the above requirements

%[row0,col0] = ind2sub(size(combined_mask),index_zeros); % convert to row column indices
[pixels.row, pixels.col] = ind2sub(size(combined_mask), pixels.index); % conver to row column indices





% check to see if we found any indexes to keep
if isempty(pixels.index)==true

    error([newline, 'No suitable pixels found!', newline])

else
    disp([newline, num2str(numel(pixels.index)),' suitable pixels found!', newline])


end


%     -----------------------------------------------------------
%     -----------------------------------------------------------



end






