% --- Create the measurement vector needed for Rodgers solution of gaussian
% model pdf and measurement pdf ---


% By Andrew J. Buggee

%%

function y = create_measurement_vector(modis, GN_inputs, pixels2use)


% We will retireve cloud properties for n pixels, designated by the
% Gauss_Newton input structure
num_pixels = GN_inputs.numPixels2Calculate;

% Use the bands specified by GN_inputs.bands2use.  Most of the time this is
% the first 7 bands

y = zeros(length(GN_inputs.bands2use),num_pixels);

for pp = 1:num_pixels
    
    
    row = pixels2use.res1km.row(pp);
    col = pixels2use.res1km.col(pp);
    
    y(:,pp) = reshape(modis.EV1km.reflectance(row,col,GN_inputs.bands2use),[],1);
    
end


end