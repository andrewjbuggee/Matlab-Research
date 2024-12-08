%% Checking retrieval uncertinaty for MODIS when using 7 MODIS channels and 35 EMIT channels
% Compare this with the King and Vaughan (2012) results from figure 7


%% Load the hyperspectral retrieval results

load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/17_Jan_2024_coast/',...
    'hyperspectral_reflectance_calculations_07-Dec-2024.mat'])

post_cov_hyperspec = retrieval.posterior_cov;



