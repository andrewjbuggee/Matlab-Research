%% Determine the cloud thermodynamic phase


% By Andrew John Buggee

%%


function [cloud_phase] = determine_cloud_phase_emit(emit, pixels2use)


% set a reflectance threshold at 500nm, below which we will consider the
% cloud to thin to decisively determine cloud phase
R_500nm_limit = 0.28;        % 1/sr



% define the number of pixels (independent spetra) being analysed
n_pix = length(pixels2use.idx);


% --- Define the spectral regions that will be used to compute percent
% differences ---

% -------------------------------------------------------------------
% ------------------ 1600 - 1700 nm Group ---------------------------
% -------------------------------------------------------------------
% find the indexes for the reflectance between 1615 and 1730 nm
idx_1700_group = find(emit.radiance.wavelength>1615 & emit.radiance.wavelength<1730);

% define the 1700 group wavelengths
wl_1700 = emit.radiance.wavelength(idx_1700_group);     % nm

% we want to store the smoothed reflectance around 1700 nm
smooth_Refl_model_1700 = zeros(length(idx_1700_group), n_pix);

% next, we use the smoothed relfectance to compute the spectral shape
% parameter, which is the percent difference of reflectance over some
% region

% find the wavelength index for the channels closest to 1.7 microns and
% 1.64 microns
[~, idx_1714] = min(abs(wl_1700 - 1714));
[~, idx_1625] = min(abs(wl_1700 - 1625));

% store the spectral shape parameter values for the 1600 micron grouping
S_1700 = zeros(length(pixels2use.idx));

% Store the reflectance at the higher end of the spectral shape parameter
R_1714 = zeros(length(pixels2use.idx));
% -------------------------------------------------------------------




% -------------------------------------------------------------------
% ------------------ 2100 - 2200 nm Group ---------------------------
% -------------------------------------------------------------------
% next store the values for the 2200 micron grouping
% find the indexes for the reflectance between 2140 and 2260 nm
idx_2200_group = find(emit.radiance.wavelength>2140 & emit.radiance.wavelength<2260);

% define the 2200 group wavelengths
wl_2200 = emit.radiance.wavelength(idx_2200_group);     % nm

% we want to store the smoothed reflectance around 1700 nm
smooth_Refl_model_2200 = zeros(length(idx_2200_group), n_pix);


% find the wavelength index for the channels closest to 2.2 microns and
% 2.15 microns
[~, idx_2240] = min(abs(wl_2200 - 2240));
[~, idx_2160] = min(abs(wl_2200 - 2160));


% store the spectral shape parameter values for the 2100 micron grouping
S_2200 = zeros(length(pixels2use.idx));

% Store the reflectance at the higher end of the spectral shape parameter
R_2240 = zeros(length(pixels2use.idx));
% -------------------------------------------------------------------





%% Step through each pixel and compute the smoothed reflectance, and spectral shape parameter

% compute the n-point moving average of the reflectance
n_avg = 4;


% step through each pixel

for pp = 1:length(pixels2use.idx)

    % Compute the smoothed reflectance over 1700 and 2200 nm
    % By smoothing, we remove small extinction features
    smooth_Refl_model_1700(:, pp) = movmean(emit.radiance.measurements(idx_1700_group, pp), n_avg);
    smooth_Refl_model_2200(:, pp) = movmean(emit.radiance.measurements(idx_2200_group, pp), n_avg);


    % store the reflectance at the higher points of the spectral shape
    % parameter equations for each spectral group
    R_1714(pp) = smooth_Refl_model_1700(idx_1714, pp);
    R_2240(pp) = smooth_Refl_model_2200(idx_2240, pp);



    % compute the spectral shape parameter (Knap et al., 2002; eq 2)
    S_1700(pp) = 100 * (smooth_Refl_model_1700(idx_1714, pp) - smooth_Refl_model_1700(idx_1625, pp))/...
        smooth_Refl_model_1700(idx_1625, pp);


    % compute the spectral shape parameter (Knap et al., 2002; eq 2)
    S_2200(pp) = 100* (smooth_Refl_model_2200(idx_2240, pp) - smooth_Refl_model_2200(idx_2160, pp))/...
        smooth_Refl_model_2200(idx_2160, pp);



end


% Sum the two shape parameters together
S = S_1700 + S_2200;


%% Use stored look-up tables to discern the thermodynamic phase

% Let's compute the L2 norm difference between the two spectral shape
% parameters and the two reflectances.

% With each



%% Let's use a very simply metric to determine if the cloud is composed of liquid droplets or solid particles

% phase identifier:
%   0 - cloud may be too thin for phase discrimination
%   1 - cloud is composed of liquid water droplets
%   2 - cloud is composed of ice particles


% if S is below 3.5, it is likely composed of liquid water droplets
% if S is greater than 3.5, it is likely composed of ice particles
liquid_solid_threshold = 3.5;

cloud_phase = zeros(1, length(pixels2use.idx));

% define the pixels that are associated with liquid water clouds
cloud_phase(S <= liquid_solid_threshold) = 1;

% define the pixels that are associated with iceclouds
cloud_phase(S > liquid_solid_threshold) = 1;

% All remaining values are undefined becase of a small reflectance at
% 500nm, indicating an optically thin cloud




end