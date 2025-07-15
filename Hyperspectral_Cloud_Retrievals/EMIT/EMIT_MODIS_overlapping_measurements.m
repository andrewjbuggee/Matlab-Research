%% Finding overlapping data between EMIT and MODIS

% By Andrew John Buggee
%%

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

results.cloudTopHeight = modis.cloud.topHeight(idx_combined_master);

results.inhomogeneity_idx = H_860(idx_combined_master);

results.pixel_time_decimal = modis.EV1km.pixel_time_decimal(idx_combined_master);
results.pixel_time_UTC = modis.EV1km.pixel_time_UTC(idx_combined_master);

% find the row and column for each
for nn = 1:numel(results.lat)

    [~, idx_min_modis] = min(abs(modis.geo.lat - results.lat(nn)), [], "all");

    [~, idx_min_emit] = min(abs(emit.radiance.geo.lat - results.lat(nn)), [], "all");

    [results.row_modis(nn), results.col_modis(nn)] = ind2sub(size(modis.cloud.aboveWaterVaporCol), idx_min_modis);

    [results.row_emit(nn), results.col_emit(nn)] = ind2sub(size(emit.radiance.measurements(:,:,1)), idx_min_emit);

end



