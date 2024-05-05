%% Remove the all pixels from the EMIT data cube that are not needed


% By Andrew John Buggee

function emit = remove_unwanted_emit_data(emit, pixels2use)


% define the number of wavelengths
num_wavelengths = size(emit.radiance.measurements, 3);

% define the number of pixels
num_pixels = length(pixels2use.row)*length(pixels2use.col);


% remove all data points except for the ones listed in pixels 2 use
keep_emit_data_index = zeros(size(emit.radiance.measurements));

% create a zero matrix for the radiance data
rad_data2keep = zeros(num_wavelengths, num_pixels);

for rr = 1:length(pixels2use.row)
    for cc = 1:length(pixels2use.col)

        keep_emit_data_index(pixels2use.row(rr), pixels2use.col(cc), :) = 1;

        % the radiance data needs to be dealt with in the for loop because
        % it is a data cube
        rad_data2keep(:, rr*cc) = reshape(emit.radiance.measurements(pixels2use.row(rr),...
            pixels2use.col(cc), :), num_wavelengths, 1);
    
    end
end

% turn the index into a logical array
keep_emit_data_index = logical(keep_emit_data_index);

% remove all unwanted data from the radiance data cube
emit.radiance.measurements = rad_data2keep;

emit.radiance.lat(~keep_emit_data_index(:,:,1)) = [];
emit.radiance.long(~keep_emit_data_index(:,:,1)) = [];

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


