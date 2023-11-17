%% ---- Estimate Cloud Depth ----



% By Andrew J. Buggee
%%

function z0 = estimate_cloud_depth(data_inputs,sensor_zenith)

% ----- parse through the inputs -----

pixels2use = data_inputs.pixels2use;
tau_c = data_inputs.truthTable.modisR17; % lets use the retrieval using bands 1 and 7.

% since band 7 is what mostly determines the optical depth, lets use this
% wavelength to determine the cloud depth

lambda = modisBandsCenter(7); % this returns the center wavelength of MODIS band 7
% perhaps in the future I can store a txt fiel with many values of the bulk
% absorption coefficient, but for now I will hard code it

k_7 = 5000; % m^(-1) - at 2130 nm the bulk absorption coefficient is about 5000 inverse meters

% -------------------------------------
% lets find all of the sensor azimuth values for the pixels that we used to
% calcualte tau

num_pixels_1km = length(pixels2use.res1km.row);
sensor_zenith_2use = zeros(num_pixels_1km,1);

for pp = 1:num_pixels_1km
    
    row = pixels2use.res1km.row(pp);
    col = pixels2use.res1km.col(pp);
    sensor_zenith_2use(pp) = sensor_zenith(row,col);
    
end

z0 = tau_c .* cosd(sensor_zenith_2use)./k_7;

end

