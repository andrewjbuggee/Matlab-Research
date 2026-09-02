%% Extinction efficiency averaged across each measured size bin
%
% DESCRIPTION:
%   Returns the Mie extinction efficiency Q_ext for each bin of a measured
%   droplet size distribution, averaged across the bin and weighted by r^2.
%
%   This is the Q that pairs with a per-bin droplet concentration in Hansen &
%   Travis (1974) Equation (2.48):
%
%       beta(z) = pi * sum_i [ Q_ext_bin(i) * r_i^2 * N_i(z) ]
%
%   Monodisperse Q_ext oscillates with a period of roughly 0.27 um in radius
%   near r = 2 um, finer than the narrowest probe bin, so evaluating it at a
%   single bin-center radius picks an arbitrary phase of that ripple: MIEV0
%   gives 2.19 at r = 6.25 um where a 1 um-resolution lookup table gives 2.05 at
%   r = 6 um. Both are correct and neither is representative. What a bin
%   actually contributes is int(Q(r) r^2 n(r) dr) across the bin, and n is only
%   known as a bin total, so the representative efficiency is the r^2-weighted
%   mean of Q across the bin.
%
%   One monodisperse Mie run covers every bin, every height and every profile,
%   so this should be called once per data set rather than per profile.
%
% INPUTS:
%   bin_edges_um    - (n_bins+1 x 1) size bin edges [microns], monotonic
%   bin_centers_um  - (n_bins x 1)   size bin centers [microns]
%   wavelength_nm   - scalar wavelength [nm]
%   which_computer  - 'anbu8374', 'andrewbuggee' or 'curc'
%
% OUTPUT:
%   Q_ext_bin       - (n_bins x 1) r^2-weighted mean extinction efficiency
%
% NOTES:
%   - Bins above r_max_um take the geometric-optics limit Q = 2. In the ORACLES
%     and VOCALS-REx profiles those bins hold less than 0.01% of the area
%     integral, so the approximation costs nothing.
%
% SEE ALSO:
%   average_mie_over_size_distribution, find_verticalProfiles_ORACLES,
%   find_verticalProfiles_VOCALS_REx_ver2
%
% By Andrew John Buggee
%%

function Q_ext_bin = binAveraged_Qext_forBins(bin_edges_um, bin_centers_um, ...
    wavelength_nm, which_computer)
% Extinction efficiency averaged across each probe size bin, weighted by r^2.
%
% Monodisperse Q_ext oscillates with a period of roughly 0.27 um in radius near
% r = 2 um, which is finer than the narrowest ORACLES bin, so evaluating it at a
% single bin-center radius picks an arbitrary phase of that ripple. What each bin
% actually contributes to Hansen & Travis (2.48) is int(Q(r) r^2 n(r) dr) across
% the bin, and n is only known as a bin total, so the representative efficiency
% is the r^2-weighted mean of Q across the bin.
%
% One monodisperse Mie run covers every bin, every height and every profile.

% Fine grid spanning the bins that carry the extinction. 200 um is well beyond
% the radius holding 99.99% of the area integral in the ORACLES profiles; the
% larger drizzle and rain bins take the geometric-optics limit.
r_min_um  = max(1.5, bin_edges_um(1));
r_max_um  = 200;
r_step_um = 0.05;
Q_LARGE_DROP = 2.0;

if strcmp(which_computer, 'anbu8374')==true
    mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

elseif strcmp(which_computer, 'andrewbuggee')==true
    mie_folder = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

elseif strcmp(which_computer, 'curc')==true
    mie_folder = '/projects/anbu8374/Mie_Calculations/';

else
    error([newline, 'Unrecognized computer name: ', which_computer, newline])
end

% --- One monodisperse Mie calculation over the fine radius grid ---
r_fine = r_min_um : r_step_um : r_max_um;              % microns

% MIEV0 rebuilds the radius vector from [start, end, step]; nudge the step until
% its vector has the same length as ours (same guard used in
% average_mie_over_size_distribution).
radius_step = mean(diff(r_fine));
while numel(r_fine(1):radius_step:r_fine(end)) ~= numel(r_fine)
    radius_step = radius_step - 1e-15;
end

[input_filename, output_filename] = write_mie_file('MIEV0', 'water', ...
    [r_fine(1), r_fine(end), radius_step], wavelength_nm, {'mono', []}, ...
    'verbose', 991, randi([0, 1000], 1), mie_folder);

[~] = runMIE(mie_folder, input_filename, output_filename, which_computer);

[ds, ~, ~] = readMIE(mie_folder, output_filename);

r_fine = r_fine(:);
Q_fine = reshape(ds.Qext(1, :), [], 1);
w_fine = r_fine.^2;

% --- Average across each probe bin ---
n_bins    = length(bin_centers_um);
Q_ext_bin = zeros(n_bins, 1);

for bb = 1:n_bins

    in_bin = r_fine >= bin_edges_um(bb) & r_fine < bin_edges_um(bb+1);

    if any(in_bin)
        Q_ext_bin(bb) = sum(Q_fine(in_bin) .* w_fine(in_bin)) / sum(w_fine(in_bin));

    elseif bin_edges_um(bb) > r_fine(end)
        % Above the Mie grid: geometric optics. These bins hold < 0.01% of the
        % area integral in every ORACLES profile checked.
        Q_ext_bin(bb) = Q_LARGE_DROP;

    else
        % Bin narrower than the grid spacing
        Q_ext_bin(bb) = interp1(r_fine, Q_fine, bin_centers_um(bb), 'linear', 'extrap');
    end

end

end
