%% Compute the FWHM of the MODIS Spectral Response Functions

% By Andrew John Buggee

%%

clear variables

% define the bands to use
bands2use = [1:7];

% define the spectral resolution
res = 0.1;  % nm

% read in the spectral response functions for the channels specified above

spec_response = modis_terra_specResponse_func(bands2use, res);

% compute the full-width at half max for each
fwhm = zeros(length(bands2use), 1);

bandwidth = zeros(length(bands2use), 2);   % nm

for nn = 1:length(bands2use)

    [max_val, idx_maxVal] = max(spec_response{nn}(:,2));

    % find the index associated with the value closest to HALF the max
    % value above and below the max value

    % below max value
    [~, min_idx_below] = min(abs(spec_response{nn}(1:idx_maxVal-1, 2) - max_val/2));

    % above max value
    [~, min_idx_above] = min(abs(spec_response{nn}(idx_maxVal+1:end, 2) - max_val/2));

    fwhm(nn) = spec_response{nn}(idx_maxVal +1 + min_idx_above,1) - spec_response{nn}(min_idx_below, 1);  % nm

    bandwidth(nn,:) = [spec_response{nn}(1,1), spec_response{nn}(end,1)];

end
