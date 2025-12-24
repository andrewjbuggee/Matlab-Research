%% Remove all pixels from the MODIS data structure that are not needed

% Pixels are defined in row and column space

% By Andrew John Buggee

%%

function modis = remove_unwanted_modis_data(modis, pixels2use)


% define the number of pixels
if isfield(pixels2use, 'idx')==true

    num_pixels = length(pixels2use.idx);

elseif isfield(pixels2use, 'linear_idx')==true

    num_pixels = length(pixels2use.linear_idx);
end


% define the rows and columns of data to keep
row = [pixels2use(:).row];
col = [pixels2use(:).col];

% Get the size of the full resolution data (2030×1354)
full_size = size(modis.geo.lat);

% Create a logical index for full resolution data
keep_modis_data_index = false(full_size);

% Mark the pixels to keep
for nn = 1:num_pixels
    keep_modis_data_index(row(nn), col(nn)) = true;
end


%% Process sensor fields (all are 2030×1354)
modis.sensor.gound_height = modis.sensor.gound_height(keep_modis_data_index);
modis.sensor.range = modis.sensor.range(keep_modis_data_index);
modis.sensor.azimuth = modis.sensor.azimuth(keep_modis_data_index);
modis.sensor.zenith = modis.sensor.zenith(keep_modis_data_index);


%% Process solar fields (all are 2030×1354)
modis.solar.azimuth = modis.solar.azimuth(keep_modis_data_index);
modis.solar.zenith = modis.solar.zenith(keep_modis_data_index);


%% Process geo fields (all are 2030×1354)
modis.geo.lat = modis.geo.lat(keep_modis_data_index);
modis.geo.long = modis.geo.long(keep_modis_data_index);


%% Process vapor fields
% col_nir and col_nir_correction are 2030×1354
modis.vapor.col_nir = modis.vapor.col_nir(keep_modis_data_index);
modis.vapor.col_nir_correction = modis.vapor.col_nir_correction(keep_modis_data_index);

% col_ir is lower resolution (406×270), needs special handling
if ~isempty(modis.vapor.col_ir)
    % Determine which lower resolution pixels correspond to the selected high-res pixels
    % Assuming 5:1 ratio (2030/406 ≈ 5, 1354/270 ≈ 5)
    low_res_size = size(modis.vapor.col_ir);
    scale_factor_row = full_size(1) / low_res_size(1);
    scale_factor_col = full_size(2) / low_res_size(2);
    
    keep_low_res_index = false(low_res_size);
    for nn = 1:num_pixels
        low_res_row = ceil(row(nn) / scale_factor_row);
        low_res_col = ceil(col(nn) / scale_factor_col);
        keep_low_res_index(low_res_row, low_res_col) = true;
    end
    
    modis.vapor.col_ir = modis.vapor.col_ir(keep_low_res_index);
end


%% Process cloud fields (all are 2030×1354)
modis.cloud.effRadius17 = modis.cloud.effRadius17(keep_modis_data_index);
modis.cloud.effRad_uncert_17 = modis.cloud.effRad_uncert_17(keep_modis_data_index);
modis.cloud.optThickness17 = modis.cloud.optThickness17(keep_modis_data_index);
modis.cloud.optThickness_uncert_17 = modis.cloud.optThickness_uncert_17(keep_modis_data_index);
modis.cloud.effRadius16 = modis.cloud.effRadius16(keep_modis_data_index);
modis.cloud.effRad_uncert_16 = modis.cloud.effRad_uncert_16(keep_modis_data_index);
modis.cloud.optThickness16 = modis.cloud.optThickness16(keep_modis_data_index);
modis.cloud.optThickness_uncert_16 = modis.cloud.optThickness_uncert_16(keep_modis_data_index);
modis.cloud.topHeight = modis.cloud.topHeight(keep_modis_data_index);
modis.cloud.topPressure = modis.cloud.topPressure(keep_modis_data_index);
modis.cloud.topTemperature = modis.cloud.topTemperature(keep_modis_data_index);
modis.cloud.phase = modis.cloud.phase(keep_modis_data_index);
modis.cloud.lwp = modis.cloud.lwp(keep_modis_data_index);
modis.cloud.aboveWaterVaporCol = modis.cloud.aboveWaterVaporCol(keep_modis_data_index);


%% Process 3D cloud fields (2030×1354×2)
% For SPI and inhomogeneity_index, keep the third dimension intact
for layer = 1:size(modis.cloud.SPI, 3)

    % the 3D structure is destroyed by the proceess. Store the original
    % data
    if layer==1
        original_SPI = modis.cloud.SPI;
        modis.cloud.SPI = [];
    end

    temp_SPI = original_SPI(:,:,layer);

    modis.cloud.SPI(:,layer) = temp_SPI(keep_modis_data_index);
end


for layer = 1:size(modis.cloud.inhomogeneity_index, 3)

    % the 3D structure is destroyed by the proceess. Store the original
    % data
    if layer==1
        original_inhomogeneity_index = modis.cloud.inhomogeneity_index;
        modis.cloud.inhomogeneity_index = [];
    end

    temp_inhom = original_inhomogeneity_index(:,:,layer);
    modis.cloud.inhomogeneity_index(:,layer) = temp_inhom(keep_modis_data_index);
end



%% Process low resolution cloud fields (406×270)
if ~isempty(modis.cloud.fraction)
    % Use the same low resolution index from vapor processing
    low_res_size = size(modis.cloud.fraction);
    scale_factor_row = full_size(1) / low_res_size(1);
    scale_factor_col = full_size(2) / low_res_size(2);
    
    keep_low_res_index = false(low_res_size);
    for nn = 1:num_pixels
        low_res_row = ceil(row(nn) / scale_factor_row);
        low_res_col = ceil(col(nn) / scale_factor_col);
        keep_low_res_index(low_res_row, low_res_col) = true;
    end
    
    modis.cloud.fraction = modis.cloud.fraction(keep_low_res_index);
    modis.cloud.fraction_day = modis.cloud.fraction_day(keep_low_res_index);
end


%% Reshape single pixel case
if num_pixels == 1
    % Reshape all vectors to column vectors for consistency
    fields_to_reshape = {'sensor', 'solar', 'geo', 'vapor', 'cloud'};
    
    for ff = 1:length(fields_to_reshape)
        field_names = fieldnames(modis.(fields_to_reshape{ff}));
        for fn = 1:length(field_names)
            data = modis.(fields_to_reshape{ff}).(field_names{fn});
            if isvector(data) && ~isempty(data)
                modis.(fields_to_reshape{ff}).(field_names{fn}) = data(:);
            end
        end
    end
end


end
