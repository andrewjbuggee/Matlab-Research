%% Remove all pixels from the EMIT data cube that are not needed

% Pixels are defined in column space

% By Andrew John Buggee

%%

function emit = remove_unwanted_emit_data(emit, pixels2use)


% define the number of wavelengths
num_wavelengths = size(emit.radiance.measurements, 3);

% define the number of pixels
if isfield(pixels2use, 'idx')==true

    num_pixels = length(pixels2use.idx);

elseif isfield(pixels2use, 'linear_idx')==true

    num_pixels = length(pixels2use.linear_idx);
end


% define the rows and columns of data to keep
row = [pixels2use(:).row];
col = [pixels2use(:).col];

% remove all data points except for the ones listed in pixels 2 run
% keep_emit_data_index = zeros(size(emit.radiance.measurements));

% create a zero matrix for the radiance data
rad_data2keep = zeros(num_wavelengths, num_pixels);

% store the geolocation data
lat_2keep = zeros(1, num_pixels);
long_2keep = zeros(1, num_pixels);

% store the solar and sensor geometry data and the time each pixel was
% recorded
sensor_az_2keep = zeros(1, num_pixels);
sensor_zen_2keep = zeros(1, num_pixels);
solar_az_2keep = zeros(1, num_pixels);
solar_zen_2keep = zeros(1, num_pixels);
utc_time_2keep = zeros(1, num_pixels);


for nn = 1:num_pixels

    % keep_emit_data_index(row(nn), col(nn), :) = 1;

    % the radiance data needs to be dealt with in the for loop because
    % it is a data cube
    rad_data2keep(:, nn) = reshape(emit.radiance.measurements(row(nn),...
        col(nn), :), num_wavelengths, 1);

    lat_2keep(:,nn) = emit.radiance.geo.lat(row(nn), col(nn));
    long_2keep(:,nn) = emit.radiance.geo.long(row(nn), col(nn));

    sensor_az_2keep(:,nn) = emit.obs.sensor.azimuth(row(nn), col(nn));
    sensor_zen_2keep(:,nn) = emit.obs.sensor.zenith(row(nn), col(nn));

    solar_az_2keep(:,nn) = emit.obs.solar.azimuth(row(nn), col(nn));
    solar_zen_2keep(:,nn) = emit.obs.solar.zenith(row(nn), col(nn));

    utc_time_2keep(:,nn) = emit.obs.utc_time(row(nn), col(nn));

end

% turn the index into a logical array
% keep_emit_data_index = logical(keep_emit_data_index);

% remove all unwanted data from the radiance data cube
emit.radiance.measurements = rad_data2keep;

emit.radiance.geo.lat = lat_2keep;
emit.radiance.geo.long = long_2keep;

% remove all unwanted data from the observation cube
emit.obs.sensor.azimuth = sensor_az_2keep;
emit.obs.sensor.zenith = sensor_zen_2keep;
emit.obs.solar.azimuth = solar_az_2keep;
emit.obs.solar.zenith = solar_zen_2keep;
emit.obs.utc_time = utc_time_2keep;

% Sort the radiance data into a desired format

if num_pixels==1

    emit.radiance.measurements = reshape(emit.radiance.measurements, [], 1);

else


end


