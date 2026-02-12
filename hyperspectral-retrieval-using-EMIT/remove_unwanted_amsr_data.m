%% Remove all pixels from the AMSR data structure that are not needed

% Pixels are defined in row and column space

% By Andrew John Buggee

%%

function AMSR = remove_unwanted_amsr_data(AMSR, pixels2use)


% define the number of pixels
if isfield(pixels2use, 'idx')==true

    num_pixels = length(pixels2use.idx);

elseif isfield(pixels2use, 'linear_idx')==true

    num_pixels = length(pixels2use.linear_idx);
end


% define the rows and columns of data to keep
row = [pixels2use(:).row];
col = [pixels2use(:).col];

% Get the size of the full resolution data (npix x nscans)
full_size = size(AMSR.geo.Latitude);

% Create a logical index for full resolution data
keep_amsr_data_index = false(full_size);

% Mark the pixels to keep
for nn = 1:num_pixels
    keep_amsr_data_index(row(nn), col(nn)) = true;
end


%% Process geolocation fields (npix x nscans)
AMSR.geo.Latitude = AMSR.geo.Latitude(keep_amsr_data_index);
AMSR.geo.Longitude = AMSR.geo.Longitude(keep_amsr_data_index);


%% Process time fields
% Time (nscans) - this is 1D per scan, need special handling
% First determine which scans we're keeping
scans_to_keep = any(keep_amsr_data_index, 2);  % Keep scan if any pixel in that scan is selected
AMSR.Time = AMSR.Time(scans_to_keep);
AMSR.Time_UTC = AMSR.Time_UTC(scans_to_keep);

% TimeHR (ntime x nscans) - keep only selected scans
AMSR.TimeHR = AMSR.TimeHR(scans_to_keep, :);


%% Process water vapor fields (npix x nscans)
AMSR.H2O.TotalPrecipitableWater = AMSR.H2O.TotalPrecipitableWater(keep_amsr_data_index);
AMSR.H2O.ErrorTPW = AMSR.H2O.ErrorTPW(keep_amsr_data_index);


%% Process cloud fields (npix x nscans)
AMSR.cloud.LiquidWaterPath = AMSR.cloud.LiquidWaterPath(keep_amsr_data_index);
AMSR.cloud.ErrorLWP = AMSR.cloud.ErrorLWP(keep_amsr_data_index);


%% Process wind fields (npix x nscans)
AMSR.wind.WindSpeed = AMSR.wind.WindSpeed(keep_amsr_data_index);
AMSR.wind.ErrorWind = AMSR.wind.ErrorWind(keep_amsr_data_index);


%% Process sea surface temperature (npix x nscans)
AMSR.SST.ReynoldsSST = AMSR.SST.ReynoldsSST(keep_amsr_data_index);


%% Process quality control fields (npix x nscans)
AMSR.QualityFlag = AMSR.QualityFlag(keep_amsr_data_index);
AMSR.ChiSquared = AMSR.ChiSquared(keep_amsr_data_index);


%% Process viewing geometry fields (if present)
% These fields may not always be present in the data structure
% if isfield(AMSR, 'solar') && isfield(AMSR.solar, 'GlintAngle')
%     AMSR.solar.GlintAngle = AMSR.solar.GlintAngle(keep_amsr_data_index);
% end


%% Process surface type fields (if present)
% if isfield(AMSR, 'LandPercentage')
%     AMSR.LandPercentage = AMSR.LandPercentage(keep_amsr_data_index);
% end


%% Process metadata fields (if they exist and are spatial)
% Note: Metadata fields are typically granule-level and don't need filtering
if isfield(AMSR, 'metadata')
    % Metadata is typically granule-level, so we keep it as-is
    % No filtering needed for scalar metadata
end


%% Reshape single pixel case
if num_pixels == 1
    % Reshape all vectors to column vectors for consistency
    
    % Geolocation
    AMSR.geo.Latitude = AMSR.geo.Latitude(:);
    AMSR.geo.Longitude = AMSR.geo.Longitude(:);
    
    % Time - these are already handled properly
    % AMSR.Time and AMSR.Time_UTC are per-scan, so they're fine as-is
    
    % Water vapor
    AMSR.H2O.TotalPrecipitableWater = AMSR.H2O.TotalPrecipitableWater(:);
    AMSR.H2O.ErrorTPW = AMSR.H2O.ErrorTPW(:);
    
    % Cloud
    AMSR.cloud.LiquidWaterPath = AMSR.cloud.LiquidWaterPath(:);
    AMSR.cloud.ErrorLWP = AMSR.cloud.ErrorLWP(:);
    
    % Wind
    AMSR.wind.WindSpeed = AMSR.wind.WindSpeed(:);
    AMSR.wind.ErrorWind = AMSR.wind.ErrorWind(:);
    
    % SST
    AMSR.SST.ReynoldsSST = AMSR.SST.ReynoldsSST(:);
    
    % Quality control
    AMSR.QualityFlag = AMSR.QualityFlag(:);
    AMSR.ChiSquared = AMSR.ChiSquared(:);
    
    % Viewing geometry (if present)
    % if isfield(AMSR, 'solar') && isfield(AMSR.solar, 'GlintAngle')
    %     AMSR.solar.GlintAngle = AMSR.solar.GlintAngle(:);
    % end
    
    % Surface type (if present)
    % if isfield(AMSR, 'LandPercentage')
    %     AMSR.LandPercentage = AMSR.LandPercentage(:);
    % end
end


end
