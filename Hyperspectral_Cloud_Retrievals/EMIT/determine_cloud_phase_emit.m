%% Determine the cloud thermodynamic phase


% By Andrew John Buggee

%%


function [cloud_phase] = determine_cloud_phase_emit(emit, pixels2use)


% define the number of pixels (independent spetra) being analysed
n_pix = length(pixels2use.idx);


% find the indexes for the reflectance between 1615 and 1730 nm
idx_1700_group = find(emit.radiance.wavelength>1615 & emit.radiance.wavelength<1730);

% define the 1700 group wavelengths
wl_1700 = emit.radiance.wavelength(idx_1700_group);     % nm

% we want to store the smoothed reflectance around 1700 nm
smooth_Refl_model_1600 = zeros(length(idx_1700_group), n_pix);

% next store the values for the 2200 micron grouping
% find the indexes for the reflectance between 2140 and 2260 nm
idx_2200_group = find(emit.radiance.wavelength>2140 & emit.radiance.wavelength<2260);

% define the 2200 group wavelengths
wl_2200 = emit.radiance.wavelength(idx_2200_group);     % nm

% we want to store the smoothed reflectance around 1700 nm
smooth_Refl_model_2200 = zeros(length(idx_2200_group), n_pix);

% find the wavelength index for the channels closest to 1.7 microns and
% 1.64 microns
[~, idx_1714] = min(abs(wl_1700 - 1700));
[~, idx_1625] = min(abs( wl_mean(inputs.idx_1600_group) - 1640));

% store the spectral shape parameter values for the 1600 micron grouping
S_1600 = zeros(length(inputs.RT.re), length(inputs.RT.tau_c));

% find the wavelength index for the channels closest to 2.2 microns and
% 2.15 microns
[~, inputs.idx_2240] = min(abs(wl_mean - 2240));
[~, inputs.idx_2160] = min(abs(wl_mean - 2160));

% need to subtract the number of spectral channels associated with the
% first group
inputs.idx_2240 = inputs.idx_2240 - length(inputs.idx_1600_group);
inputs.idx_2160 = inputs.idx_2160 - length(inputs.idx_1600_group);

% store the spectral shape parameter values for the 2100 micron grouping
S_2100 = zeros(length(inputs.RT.re), length(inputs.RT.tau_c));



% step through each pixel

for pp = 1:length(pixels2use.idx)








% compute the n-point moving average of the reflectance
n = 4;

smooth_Refl_model_1600(rr, tc, :) = movmean(Refl_model(rr, tc, inputs.idx_1600_group), 4);
smooth_Refl_model_2100(rr, tc, :) = movmean(Refl_model(rr, tc, inputs.idx_2100_group), 4);




% compute the spectral shape parameter (Knap et al., 2002; eq 2)
S_1600(rr, tc) = 100* (smooth_Refl_model_1600(rr, tc, inputs.idx_1714) - smooth_Refl_model_1600(rr, tc, inputs.idx_1625))/...
    smooth_Refl_model_1600(rr, tc, inputs.idx_1625);


% compute the spectral shape parameter (Knap et al., 2002; eq 2)
S_2100(rr, tc) = 100* (smooth_Refl_model_2100(rr, tc, inputs.idx_2240) - smooth_Refl_model_2100(rr, tc, inputs.idx_2160))/...
    smooth_Refl_model_2100(rr, tc, inputs.idx_2160);



end

