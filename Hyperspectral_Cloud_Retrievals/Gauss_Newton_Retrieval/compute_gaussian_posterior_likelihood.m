%% --- Compute the Posterior Covariance Matrix ---


% By Andrew J. Buggee
%%

function [S_posterior, likelihood] = compute_gaussian_posterior_likelihood(inputs,data_inputs,K,modis)

% --- Read Bayesian inputs ---
model_covaraince = inputs.model.covaraince;
measurement_covaraince = inputs.measurement.covaraince;

% set up measurements into a vector
% --** For pixel 1 first **--
pixels2use = data_inputs.pixels2use;
bands2run = data_inputs.inputs.bands2run;

numPixels_500 = size(pixels2use.res500m.row,1); % number of pixels to compute in the 500 meter resolution data
num_bands = size(K,1); % number of bands
measured_data = cat(3,modis.EV.m250.reflectance,modis.EV.m500.reflectance); % combine data so all bands are easy to grab and manipulate
measured_data = measured_data(:,:,bands2run); % only use the bands that are available in the lookup table

measurements = zeros(num_bands,numPixels_500);

for pp = 1:size(pixels2use.res500m.row,1)
    
    row = pixels2use.res500m.row(pp);
    col = pixels2use.res500m.col(pp);
    
    measurements(:,pp) = reshape(measured_data(row,col,:),num_bands,1);
    
    
    
    
end







end