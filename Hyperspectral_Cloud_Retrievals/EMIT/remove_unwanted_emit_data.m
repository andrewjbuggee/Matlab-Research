%% Remove all pixels from the EMIT data cube that are not needed

% Pixels are defined in column space

% By Andrew John Buggee

%%

function emit = remove_unwanted_emit_data(emit, pixels2use)


% define the number of wavelengths
num_wavelengths = size(emit.radiance.measurements, 3);

% define the number of pixels
num_pixels = length(pixels2use.idx);

% define the rows and columns of data to keep
row = pixels2use.row;
col = pixels2use.col;

% remove all data points except for the ones listed in pixels 2 run
keep_emit_data_index = zeros(size(emit.radiance.measurements));

% create a zero matrix for the radiance data
rad_data2keep = zeros(num_wavelengths, num_pixels);

for nn = 1:num_pixels

    keep_emit_data_index(row(nn), col(nn), :) = 1;

    % the radiance data needs to be dealt with in the for loop because
    % it is a data cube
    rad_data2keep(:, nn) = reshape(emit.radiance.measurements(row(nn),...
        col(nn), :), num_wavelengths, 1);


end

% turn the index into a logical array
keep_emit_data_index = logical(keep_emit_data_index);

% remove all unwanted data from the radiance data cube
emit.radiance.measurements = rad_data2keep;

emit.radiance.geo.lat(~keep_emit_data_index(:,:,1)) = [];
emit.radiance.geo.long(~keep_emit_data_index(:,:,1)) = [];

% remove all unwanted data from the observation cube
emit.obs.sensor.azimuth(~keep_emit_data_index(:,:,1)) = [];
emit.obs.sensor.zenith(~keep_emit_data_index(:,:,1)) = [];
emit.obs.solar.azimuth(~keep_emit_data_index(:,:,1)) = [];
emit.obs.solar.zenith(~keep_emit_data_index(:,:,1)) = [];
emit.obs.utc_time(~keep_emit_data_index(:,:,1)) = [];

% Sort the radiance data into a desired format

if length(pixels2use.idx)==1

    emit.radiance.measurements = reshape(emit.radiance.measurements, [], 1);

else

end


