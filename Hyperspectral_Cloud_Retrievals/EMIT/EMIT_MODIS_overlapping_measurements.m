%% Finding overlapping data between EMIT and MODIS

% By Andrew John Buggee

%% Load EMIT DATA

clear variables

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

% emitDataFolder = '17_Jan_2024_coast/';

% 27 january has overlap with MODIS observations
emitDataFolder = '27_Jan_2024/';


% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders();

% Load EMIT Data

[emit,L1B_fileName] = retrieveEMIT_data([folder_paths.emitDataPath, emitDataFolder]);

%% Load MODIS Data

% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ***** Define the MODIS Folder *****

    modisFolder = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


elseif strcmp(whatComputer,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % ----- Define the MODIS folder name -----

    modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/MODIS_data/'];


elseif strcmp(whatComputer,'curc')==true



    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


    % Define the MODIS folder name

    modisFolder = '/projects/anbu8374/MODIS_data/';


end


% ----- January 27th 2024 at UTC 1500 -----
modisData = '2024_01_27/';


[modis,L1B_fileName] = retrieveMODIS_data([modisFolder, modisData]);


%% Grab the MODIS and EMIT lat

% Plot MODIS effective radius data using geoscatter

modis_lat = double(modis.geo.lat(:));
modis_long = double(modis.geo.long(:));
[isOcean] = land_or_ocean(modis_lat, modis_long, 10, true);
figure; 
sm = geoscatter(modis_lat(isOcean), modis_long(isOcean), 100, modis.cloud.effRadius17(isOcean), '.');
hold on
% add EMIT reflectance data on top!

emit_lat = emit.radiance.geo.lat(:);
emit_long = emit.radiance.geo.long(:);
se = geoscatter(emit_lat, emit_long, 100, 'k', 'filled');

% set the transparency low
se.MarkerFaceAlpha = 0.0055;
se.MarkerEdgeAlpha = 0.0055;


%%

% Plot MODIS inhomogeneity index data using geoscatter

modis_lat = double(modis.geo.lat(:));
modis_long = double(modis.geo.long(:));
H_860 = modis.cloud.inhomogeneity_index(:,:,2);

emit_lat = emit.radiance.geo.lat(:);
emit_long = emit.radiance.geo.long(:);

[isOcean] = land_or_ocean(modis_lat, modis_long, 10, true);


figure; 
sm = geoscatter(modis_lat(isOcean), modis_long(isOcean), 100, H_860(isOcean), '.');
hold on
colorbar
% add EMIT reflectance data on top!
se = geoscatter(emit_lat, emit_long, 100, 'k', 'filled');

% set the transparency low
se.MarkerFaceAlpha = 0.0055;
se.MarkerEdgeAlpha = 0.0055;


%% Grab MODIS data within the EMIT data set


emit_lat_min = min(emit_lat);
emit_lat_max = max(emit_lat);

emit_long_min = min(emit_long);
emit_long_max = max(emit_long);

idx_lat = modis.geo.lat>=emit_lat_min & modis.geo.lat<=emit_lat_max;
idx_long = modis.geo.long>=emit_long_min & modis.geo.long<=emit_long_max;
idx_combined = logical(idx_lat .* idx_long);

% grab only the combined values
modis_lat_emitOverlap = modis.geo.lat(idx_combined);
modis_long_emitOverlap = modis.geo.long(idx_combined);

% grab the effective radius values
H_860_emitOverlap = H_860(idx_combined);



% plot just the overlapping data
figure; 
sm2 = geoscatter(modis_lat_emitOverlap, modis_long_emitOverlap, 100, H_860_emitOverlap, '.');
hold on
colorbar

%%

% Find all MODIS pixels within the footprint of EMIT (** not exact **)
emit_lat_min = min(emit.radiance.geo.lat, [], "all");
emit_lat_max = max(emit.radiance.geo.lat, [], "all");

emit_long_min = min(emit.radiance.geo.long, [], "all");
emit_long_max = max(emit.radiance.geo.long, [], "all");

idx_lat = modis.geo.lat>=emit_lat_min & modis.geo.lat<=emit_lat_max;
idx_long = modis.geo.long>=emit_long_min & modis.geo.long<=emit_long_max;

%% Find distance between each MODIS and EMIT pixel

% we will be computing the arclength between points on an ellipsoid
% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
% wgs84 = wgs84Ellipsoid("m");
% 
% modis_lat = double(modis.geo.lat(:));
% modis_long = double(modis.geo.long(:));
% 
% emit_lat = emit.radiance.geo.lat(:);
% emit_long = emit.radiance.geo.long(:);

%% Set up an empty array for each vocals-rex profile
% modis_minDist = zeros(length(emit_long), 1);
% 
% 
% parfor nn = 1:numel(emit_long)
% 
%     dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, emit_lat(nn),emit_long(nn), wgs84);
% 
%     [modis_minDist(nn), index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');            % m - minimum distance
% 
% end




%%


% find the subset where a water cloud is present
is_liquidWater = modis.cloud.phase==2;           % 2 is the value designated for liquid water

% find a cloud with an optical thickness of atleast 4
is_tauGreaterThan4 = modis.cloud.optThickness17>=4;

% find a pixel where the effective radius retreival uncertainty is less
% than 10 percent
is_reUncertLessThan10Percent = modis.cloud.effRad_uncert_17<=10;

% find where the inhomogeneity index is less than 0.3
is_H_lessThan = modis.cloud.inhomogeneity_index(:,:,2)<1;

% let's combine these to find all values that meet these requirements
idx_combined_master = logical(idx_lat .* idx_long .* is_liquidWater .* is_tauGreaterThan4 .*...
                        is_reUncertLessThan10Percent .* is_H_lessThan);


%%




% print the results

results.lat = double(modis.geo.lat(idx_combined_master));
results.long = double(modis.geo.long(idx_combined_master));

results.effRad17 = modis.cloud.effRadius17(idx_combined_master);
results.effRad_uncert_17 = modis.cloud.effRad_uncert_17(idx_combined_master);

results.optThickness17 = modis.cloud.optThickness17(idx_combined_master);
results.optThickness_uncert_17 = modis.cloud.optThickness_uncert_17(idx_combined_master);

results.aboveWaterVaporCol = modis.cloud.aboveWaterVaporCol(idx_combined_master);

results.aboveWaterVaporCol_NIR = modis.vapor.col_nir(idx_combined_master);
results.aboveWaterVaporCol_NIR_correction = modis.vapor.col_nir_correction(idx_combined_master);

results.cloudTopHeight = modis.cloud.topHeight(idx_combined_master);

results.inhomogeneity_idx = modis.cloud.inhomogeneity_index(logical(cat(3,zeros(size(idx_combined_master)), idx_combined_master)));

results.pixel_time_decimal = modis.EV1km.pixel_time_decimal(idx_combined_master);
results.pixel_time_UTC = modis.EV1km.pixel_time_UTC(idx_combined_master);

% find the row and column for each
for nn = 1:numel(results.lat)

    % dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, results.lat(nn), results.long(nn), wgs84);
    % 
    % [modis_minDist(nn), index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');            % m - minimum distance

    [~, idx_min_modis] = min(sqrt((modis.geo.lat - results.lat(nn)).^2 + ...
        (modis.geo.long - results.long(nn)).^2), [], "all");

    [~, idx_min_emit] = min(sqrt((emit.radiance.geo.lat - results.lat(nn)).^2 + ...
        (emit.radiance.geo.long - results.long(nn)).^2), [], "all");


    [results.row_modis(nn), results.col_modis(nn)] = ind2sub(size(modis.cloud.aboveWaterVaporCol), idx_min_modis);

    [results.row_emit(nn), results.col_emit(nn)] = ind2sub(size(emit.radiance.measurements(:,:,1)), idx_min_emit);

end



