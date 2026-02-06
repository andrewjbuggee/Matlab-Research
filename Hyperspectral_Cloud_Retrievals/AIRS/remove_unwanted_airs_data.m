%% Remove all pixels from the AIRS data structure that are not needed

% Pixels are defined in row and column space

% By Andrew John Buggee

%%

function airs = remove_unwanted_airs_data(airs, pixels2use)


% define the number of pixels
if isfield(pixels2use, 'idx')==true

    num_pixels = length(pixels2use.idx);

elseif isfield(pixels2use, 'linear_idx')==true

    num_pixels = length(pixels2use.linear_idx);
end


% define the rows and columns of data to keep
row = [pixels2use(:).row];
col = [pixels2use(:).col];

% Get the size of the full resolution data (GeoTrack x GeoXTrack)
full_size = size(airs.geo.Latitude);

% Create a logical index for full resolution data
keep_airs_data_index = false(full_size);

% Mark the pixels to keep
for nn = 1:num_pixels
    keep_airs_data_index(row(nn), col(nn)) = true;
end


%% Process geolocation fields (GeoTrack x GeoXTrack)
airs.geo.Latitude = airs.geo.Latitude(keep_airs_data_index);
airs.geo.Longitude = airs.geo.Longitude(keep_airs_data_index);


%% Process time fields (GeoTrack x GeoXTrack)
airs.Time = airs.Time(keep_airs_data_index);
airs.Time_UTC = airs.Time_UTC(keep_airs_data_index);


%% Process temperature fields
% Temperature profiles (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.temp.prof_std, 3)
    
    % Store the original 3D structure
    if layer == 1
        original_temp_prof = airs.temp.prof_std;
        original_temp_prof_err = airs.temp.prof_std_err;
        airs.temp.prof_std = [];
        airs.temp.prof_std_err = [];
    end
    
    temp_slice = original_temp_prof(:,:,layer);
    airs.temp.prof_std(:,layer) = temp_slice(keep_airs_data_index);
    
    temp_err_slice = original_temp_prof_err(:,:,layer);
    airs.temp.prof_std_err(:,layer) = temp_err_slice(keep_airs_data_index);
end

% Temperature QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.temp.prof_std_QC, 3)
    
    if layer == 1
        original_temp_QC = airs.temp.prof_std_QC;
        airs.temp.prof_std_QC = [];
    end
    
    temp_QC_slice = original_temp_QC(:,:,layer);
    airs.temp.prof_std_QC(:,layer) = temp_QC_slice(keep_airs_data_index);
end

% Surface air temperature (GeoTrack x GeoXTrack)
airs.temp.SurfAir = airs.temp.SurfAir(keep_airs_data_index);
airs.temp.SurfAir_QC = airs.temp.SurfAir_QC(keep_airs_data_index);


%% Process pressure fields
% Surface pressure (GeoTrack x GeoXTrack)
airs.pressure.SurfStd = airs.pressure.SurfStd(keep_airs_data_index);
airs.pressure.SurfStd_QC = airs.pressure.SurfStd_QC(keep_airs_data_index);

% Note: pressStd and pressH2O are 1D arrays (pressure levels), not spatial, so don't filter them


%% Process Geopotential Height fields
% Temperature profiles (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.Z.GP_Height, 3)
    
    % Store the original 3D structure
    if layer == 1
        original_Z_prof = airs.Z.GP_Height;
        airs.Z.GP_Height = [];
    end
    
    Z_slice = original_Z_prof(:,:,layer);
    airs.Z.GP_Height(:,layer) = Z_slice(keep_airs_data_index);
    
end

% Z QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.Z.GP_Height_QC, 3)
    
    if layer == 1
        original_Z_QC = airs.Z.GP_Height_QC;
        airs.Z.GP_Height_QC = [];
    end
    
    Z_QC_slice = original_Z_QC(:,:,layer);
    airs.Z.GP_Height_QC(:,layer) = Z_QC_slice(keep_airs_data_index);
end

% Surface geopotential height (GeoTrack x GeoXTrack)
airs.Z.GP_Surface = airs.Z.GP_Surface(keep_airs_data_index);
airs.Z.GP_Surface_QC = airs.Z.GP_Surface_QC(keep_airs_data_index);


%% Process water vapor fields
% Water vapor mass mixing ratio (GeoTrack x GeoXTrack x H2OPressureLay)
for layer = 1:size(airs.H2O.MMR_Std, 3)
    
    if layer == 1
        original_H2O_MMR = airs.H2O.MMR_Std;
        original_H2O_MMR_err = airs.H2O.MMR_StdErr;
        airs.H2O.MMR_Std = [];
        airs.H2O.MMR_StdErr = [];
    end
    
    H2O_slice = original_H2O_MMR(:,:,layer);
    airs.H2O.MMR_Std(:,layer) = H2O_slice(keep_airs_data_index);
    
    H2O_err_slice = original_H2O_MMR_err(:,:,layer);
    airs.H2O.MMR_StdErr(:,layer) = H2O_err_slice(keep_airs_data_index);
end

% Water vapor QC (GeoTrack x GeoXTrack x H2OPressureLay)
for layer = 1:size(airs.H2O.MMR_Std_QC, 3)
    
    if layer == 1
        original_H2O_QC = airs.H2O.MMR_Std_QC;
        airs.H2O.MMR_Std_QC = [];
    end
    
    H2O_QC_slice = original_H2O_QC(:,:,layer);
    airs.H2O.MMR_Std_QC(:,layer) = H2O_QC_slice(keep_airs_data_index);
end

% Water vapor at pressure levels (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.H2O.MMR_LevStd, 3)
    
    if layer == 1
        original_H2O_Lev = airs.H2O.MMR_LevStd;
        airs.H2O.MMR_LevStd = [];
    end
    
    H2O_Lev_slice = original_H2O_Lev(:,:,layer);
    airs.H2O.MMR_LevStd(:,layer) = H2O_Lev_slice(keep_airs_data_index);
end

% Water vapor level QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.H2O.MMR_LevStd_QC, 3)
    
    if layer == 1
        original_H2O_Lev_QC = airs.H2O.MMR_LevStd_QC;
        airs.H2O.MMR_LevStd_QC = [];
    end
    
    H2O_Lev_QC_slice = original_H2O_Lev_QC(:,:,layer);
    airs.H2O.MMR_LevStd_QC(:,layer) = H2O_Lev_QC_slice(keep_airs_data_index);
end

% Total column water vapor (GeoTrack x GeoXTrack)
airs.H2O.totCol_Std = airs.H2O.totCol_Std(keep_airs_data_index);
airs.H2O.totCol_Std_QC = airs.H2O.totCol_Std_QC(keep_airs_data_index);

% Relative humidity (GeoTrack x GeoXTrack x H2OPressureLay)
for layer = 1:size(airs.H2O.RelHum, 3)
    
    if layer == 1
        original_RelHum = airs.H2O.RelHum;
        airs.H2O.RelHum = [];
    end
    
    RelHum_slice = original_RelHum(:,:,layer);
    airs.H2O.RelHum(:,layer) = RelHum_slice(keep_airs_data_index);
end

% Relative humidity QC (GeoTrack x GeoXTrack x H2OPressureLay)
for layer = 1:size(airs.H2O.RelHum_QC, 3)
    
    if layer == 1
        original_RelHum_QC = airs.H2O.RelHum_QC;
        airs.H2O.RelHum_QC = [];
    end
    
    RelHum_QC_slice = original_RelHum_QC(:,:,layer);
    airs.H2O.RelHum_QC(:,layer) = RelHum_QC_slice(keep_airs_data_index);
end


%% Process methane fields
% Methane VMR profiles (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.CH4.VMR_LevStd, 3)
    
    if layer == 1
        original_CH4_VMR = airs.CH4.VMR_LevStd;
        original_CH4_VMR_err = airs.CH4.VMR_LevStdErr;
        airs.CH4.VMR_LevStd = [];
        airs.CH4.VMR_LevStdErr = [];
    end
    
    CH4_slice = original_CH4_VMR(:,:,layer);
    airs.CH4.VMR_LevStd(:,layer) = CH4_slice(keep_airs_data_index);
    
    CH4_err_slice = original_CH4_VMR_err(:,:,layer);
    airs.CH4.VMR_LevStdErr(:,layer) = CH4_err_slice(keep_airs_data_index);
end

% Methane QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.CH4.VMR_LevStd_QC, 3)
    
    if layer == 1
        original_CH4_QC = airs.CH4.VMR_LevStd_QC;
        airs.CH4.VMR_LevStd_QC = [];
    end
    
    CH4_QC_slice = original_CH4_QC(:,:,layer);
    airs.CH4.VMR_LevStd_QC(:,layer) = CH4_QC_slice(keep_airs_data_index);
end

% Total column methane (GeoTrack x GeoXTrack)
airs.CH4.total_column = airs.CH4.total_column(keep_airs_data_index);
airs.CH4.total_column_QC = airs.CH4.total_column_QC(keep_airs_data_index);

% Degrees of freedom for methane (GeoTrack x GeoXTrack)
airs.CH4.dof = airs.CH4.dof(keep_airs_data_index);


%% Process carbon monoxide fields
% CO VMR profiles (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.CO.VMR_LevStd, 3)
    
    if layer == 1
        original_CO_VMR = airs.CO.VMR_LevStd;
        original_CO_VMR_err = airs.CO.VMR_LevStdErr;
        airs.CO.VMR_LevStd = [];
        airs.CO.VMR_LevStdErr = [];
    end
    
    CO_slice = original_CO_VMR(:,:,layer);
    airs.CO.VMR_LevStd(:,layer) = CO_slice(keep_airs_data_index);
    
    CO_err_slice = original_CO_VMR_err(:,:,layer);
    airs.CO.VMR_LevStdErr(:,layer) = CO_err_slice(keep_airs_data_index);
end

% CO QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.CO.VMR_LevStd_QC, 3)
    
    if layer == 1
        original_CO_QC = airs.CO.VMR_LevStd_QC;
        airs.CO.VMR_LevStd_QC = [];
    end
    
    CO_QC_slice = original_CO_QC(:,:,layer);
    airs.CO.VMR_LevStd_QC(:,layer) = CO_QC_slice(keep_airs_data_index);
end

% Total column CO (GeoTrack x GeoXTrack)
airs.CO.total_column = airs.CO.total_column(keep_airs_data_index);
airs.CO.total_column_QC = airs.CO.total_column_QC(keep_airs_data_index);


%% Process ozone fields
% Ozone VMR profiles (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.O3.VMR_Std, 3)
    
    if layer == 1
        original_O3_VMR = airs.O3.VMR_Std;
        airs.O3.VMR_Std = [];
    end
    
    O3_slice = original_O3_VMR(:,:,layer);
    airs.O3.VMR_Std(:,layer) = O3_slice(keep_airs_data_index);
end

% Ozone QC (GeoTrack x GeoXTrack x StdPressureLev)
for layer = 1:size(airs.O3.VMR_Std_QC, 3)
    
    if layer == 1
        original_O3_QC = airs.O3.VMR_Std_QC;
        airs.O3.VMR_Std_QC = [];
    end
    
    O3_QC_slice = original_O3_QC(:,:,layer);
    airs.O3.VMR_Std_QC(:,layer) = O3_QC_slice(keep_airs_data_index);
end

% Total column ozone (GeoTrack x GeoXTrack)
airs.O3.totCol_Std = airs.O3.totCol_Std(keep_airs_data_index);
airs.O3.totCol_Std_QC = airs.O3.totCol_Std_QC(keep_airs_data_index);


%% Process cloud fields (GeoTrack x GeoXTrack)
airs.cloud.FracTot = airs.cloud.FracTot(keep_airs_data_index);
airs.cloud.FracTot_QC = airs.cloud.FracTot_QC(keep_airs_data_index);
airs.cloud.presTop_ = airs.cloud.presTop_(keep_airs_data_index);
airs.cloud.presTop_QC = airs.cloud.presTop_QC(keep_airs_data_index);


%% Process viewing geometry fields (GeoTrack x GeoXTrack)
airs.sensor.zen = airs.sensor.zen(keep_airs_data_index);
airs.solar.zen = airs.solar.zen(keep_airs_data_index);


%% Process land fraction (GeoTrack x GeoXTrack)
airs.landFrac = airs.landFrac(keep_airs_data_index);


%% Process metadata fields (if they exist and are spatial)
% Note: Most metadata fields are scalar/granule-level and don't need filtering
% Only filter if they have spatial dimensions
if isfield(airs, 'metadata')
    % Metadata is typically granule-level, so we keep it as-is
    % No filtering needed for scalar metadata
end


%% Reshape single pixel case
if num_pixels == 1
    % Reshape all vectors to column vectors for consistency
    
    % Geolocation
    airs.geo.Latitude = airs.geo.Latitude(:);
    airs.geo.Longitude = airs.geo.Longitude(:);
    
    % Time
    airs.Time = airs.Time(:);
    % Time_UTC is datetime, handled differently
    
    % Temperature
    airs.temp.SurfAir = airs.temp.SurfAir(:);
    airs.temp.SurfAir_QC = airs.temp.SurfAir_QC(:);
    
    % Pressure
    airs.pressure.SurfStd = airs.pressure.SurfStd(:);
    airs.pressure.SurfStd_QC = airs.pressure.SurfStd_QC(:);
    
    % Water vapor scalars
    airs.H2O.totCol_Std = airs.H2O.totCol_Std(:);
    airs.H2O.totCol_Std_QC = airs.H2O.totCol_Std_QC(:);
    
    % Methane scalars
    airs.CH4.total_column = airs.CH4.total_column(:);
    airs.CH4.total_column_QC = airs.CH4.total_column_QC(:);
    airs.CH4.dof = airs.CH4.dof(:);
    
    % CO scalars
    airs.CO.total_column = airs.CO.total_column(:);
    airs.CO.total_column_QC = airs.CO.total_column_QC(:);
    
    % Ozone scalars
    airs.O3.totCol_Std = airs.O3.totCol_Std(:);
    airs.O3.totCol_Std_QC = airs.O3.totCol_Std_QC(:);
    
    % Cloud
    airs.cloud.FracTot = airs.cloud.FracTot(:);
    airs.cloud.FracTot_QC = airs.cloud.FracTot_QC(:);
    airs.cloud.presTop_ = airs.cloud.presTop_(:);
    airs.cloud.presTop_QC = airs.cloud.presTop_QC(:);
    
    % Geometry
    airs.sensor.zen = airs.sensor.zen(:);
    airs.solar.zen = airs.solar.zen(:);
    
    % Land fraction
    airs.landFrac = airs.landFrac(:);
end


end
