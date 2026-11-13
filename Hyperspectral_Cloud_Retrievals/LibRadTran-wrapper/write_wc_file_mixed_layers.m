%% Write a libRadtran water-cloud (wc) file whose cloud spans exactly z_base ... z_top
%
% Drop-in replacement for write_wc_file_from_in_situ.m (same 10 inputs, same
% output) for profiles sampled at discrete levels, e.g. the 7-level synthetic
% profiles of hysics_refl_from_synthetic_NN_inputs.m.
%
% WHY A NEW WRITER
% libRadtran (>= 1.4) reads every wc_file line as a LAYER quantity: the line
% at altitude z_i sets constant LWC and r_e from z_i UP TO the next line above
% (libRadtran 2.0.6 user guide, wc_file: "interpreted as layer quantities").
% write_wc_file_from_in_situ.m writes the N sampled levels as N lines and
% closes the cloud with a zero line at z_top + dz. The top sample therefore
% fills an extra layer [z_top, z_top + dz] ABOVE the profile, and libRadtran's
% cloud is one level spacing thicker than the profile the tau_c and LWP
% labels are integrated over (trapezoid, z_base ... z_top). On 30 synthetic
% clouds this raised libRadtran's tau by a median of 7.8% (5th-95th percentile
% 1.4-24%) and its LWP by 8.2% (1.5-27%)
% (Testing/test_wc_geometry_optical_depth.m, 11 Nov 2026).
%
% WHAT THIS WRITER DOES
% N sampled levels bound N-1 intervals. Each interval [z_{k+1}, z_k] becomes
% ONE homogeneous slab holding the well-mixed contents of its two bounding
% levels:
%       LWC_mix = (LWC_k + LWC_{k+1}) / 2
%       r_e,mix = LWC_mix / [ (LWC_k / r_e,k + LWC_{k+1} / r_e,k+1) / 2 ]
% The slab's line is written at its bottom, z_{k+1}, and a zero line at
% z_top closes the cloud, so libRadtran's cloud spans exactly z_base ... z_top.
%
% SOURCE OF THE MIXING RULE: this rule is Claude's own derivation (Nov 2026),
% not a result stated in Hansen & Travis (1974). Only the definitions it uses
% are theirs: r_e = <r^3>/<r^2> (H&T eq. 2.53) and beta = int pi r^2 Q_ext n dr
% (H&T eq. 2.48). The derivation: mixing two equal volumes averages their size
% distributions, so the moments average; LWC is proportional to <r^3> and
% LWC/r_e to <r^2> (the droplet cross-section), giving the r_e,mix above.
% Full write-up, worked example and libRadtran check: Claude Doc
% "Mixed-layer r_e for libRadtran wc files"
% (https://claude.ai/artifact/AHaW9FGo7NAnvouBLmvbxT).
%
% WHAT THIS GUARANTEES
%   - LWP: libRadtran's column water equals the trapezoid LWP of the sampled
%     profile exactly (up to the printed precision).
%   - tau: equals the trapezoid of beta(z) up to the ratio
%     Qbar(r_e,mix) / [cross-section-weighted mean of Qbar_k, Qbar_k+1],
%     because libRadtran represents the mixture by one gamma distribution.
%     On 30 synthetic clouds at 550 nm: within +/-0.4%.
%   - The slab r_e values are mixtures, NOT the sampled values. The
%     radiances therefore come from a slightly smoothed version of the
%     sampled r_e profile.
%
% INPUTS (identical to write_wc_file_from_in_situ.m):
%   (1) re  - effective radius at each sampled level [microns]
%   (2) lwc - liquid water content at each sampled level [g m^-3]
%   (3) z   - altitude of each sampled level [km]; either order (cloud top
%             first or cloud base first). Levels must be strictly monotonic.
%   (4) inSitu_campaign - string used in the filename
%   (5) inSitu_date     - date string used in the filename
%   (6) inSitu_time     - number used in the filename
%   (7) compute_weighting_functions - must be false. The Platnick (2000)
%       layer-by-layer files are not implemented for this geometry.
%   (8) computer_name   - unused; kept so the call matches the old writer
%   (9) index           - unique integer tag for the filename
%   (10) wc_folder_path - folder for the .DAT file (with trailing separator)
%
% OUTPUT:
%   fileName - name (not path) of the .DAT file written in wc_folder_path
%
% By Andrew John Buggee and Claude, 11 Nov 2026
%%

function [fileName] = write_wc_file_mixed_layers(re, lwc, z, ...
    inSitu_campaign, inSitu_date, inSitu_time, compute_weighting_functions, ...
    computer_name, index, wc_folder_path) %#ok<INUSD>  computer_name unused

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

if nargin ~= 10
    error(['Should be 10 inputs (same as write_wc_file_from_in_situ): re, lwc, z, ', ...
        'campaign name, date, time, compute_weighting_functions flag, computer name, ', ...
        'file index, and the wc folder path.'])
end

if compute_weighting_functions
    error(['write_wc_file_mixed_layers does not implement the weighting-function ', ...
        'files (one file per removed bottom layer). Use write_wc_file_from_in_situ ', ...
        'for that, or extend this function.'])
end

% column vectors throughout
re  = double(re(:));    % microns
lwc = double(lwc(:));   % g m^-3
z   = double(z(:));     % km

if numel(re) ~= numel(lwc) || numel(re) ~= numel(z)
    error('The inputs re, lwc and z must be the same length.')
end
if numel(z) < 2
    error('At least two sampled levels are needed to define one cloud slab.')
end
if any(isnan(re)) || any(isnan(lwc)) || any(isnan(z))
    error('re, lwc and z must not contain NaN values.')
end
if any(lwc < 0)
    error('LWC must be non-negative.')
end

% Order the levels cloud top -> cloud base (z decreasing), whatever order
% they arrive in (an ascending or a descending aircraft profile)
[z, order] = sort(z, 'descend');
re  = re(order);
lwc = lwc(order);
if any(diff(z) == 0)
    error('Two sampled levels share the same altitude; slabs would have zero thickness.')
end

% ------------------------------------------------------------
% ------------- MIX ADJACENT LEVELS INTO SLABS ---------------
% ------------------------------------------------------------
% Slab k spans [z(k+1), z(k)], k = 1 ... N-1 (k = 1 is the top slab).

% Liquid water adds as mass, so the slab's LWC is the mean of its two levels.
% LWC_mix * thickness is then exactly the trapezoid term for that interval.
lwc_mix = (lwc(1:end-1) + lwc(2:end)) / 2;                      % g m^-3, (N-1, 1)

% Droplet cross-section adds too. LWC / r_e = (4/3) rho_w G, with G the
% geometric cross-section per unit volume (from H&T eq. 2.53), so the slab's
% cross-section factor is the mean of the two levels' LWC / r_e.
% A level with r_e = 0 has no droplets (LWC must then be 0): it contributes no
% cross-section, so 0/0 is replaced by 0.
lwc_over_re          = lwc ./ re;                               % g m^-3 um^-1
lwc_over_re(re <= 0) = 0;
xsec_mix = (lwc_over_re(1:end-1) + lwc_over_re(2:end)) / 2;     % g m^-3 um^-1

% r_e of the mixture = total <r^3> / total <r^2>  =  LWC_mix / xsec_mix
re_mix = zeros(size(lwc_mix));                                  % um
has_cloud = lwc_mix > 0 & xsec_mix > 0;
re_mix(has_cloud) = lwc_mix(has_cloud) ./ xsec_mix(has_cloud);
lwc_mix(~has_cloud) = 0;                                        % an empty slab stays empty

% ------------------------------------------------------------
% ------------------- LINES OF THE FILE ----------------------
% ------------------------------------------------------------
% libRadtran: each line holds from its altitude UP TO the next line above.
%   z(1) = z_top : zeros, closes the cloud at the profile's top
%   z(2) ... z(N): slab k = 1 ... N-1, each written at its BOTTOM, z(k+1)
%   0 km         : zeros, clear air from the surface up to cloud base
%                  (omitted when the cloud base is the surface itself)
z_lines   = [z(1);  z(2:end)];
lwc_lines = [0;     lwc_mix];
re_lines  = [0;     re_mix];
if z(end) > 0
    z_lines   = [z_lines;   0];
    lwc_lines = [lwc_lines; 0];
    re_lines  = [re_lines;  0];
end

% ------------------------------------------------------------
% ---------------------- WRITE WC FILE -----------------------
% ------------------------------------------------------------

fileName = ['WC_mixedLayers_', inSitu_campaign, '_', char(inSitu_date), '_', ...
    num2str(inSitu_time), '-UTC_nn', num2str(index), '.DAT'];

fileID = fopen([wc_folder_path, fileName], 'w');
if fileID < 0
    error('Could not open %s for writing.', [wc_folder_path, fileName])
end

% Header lines, as in write_wc_file_from_in_situ.m
fprintf(fileID, '%s %10s %9s %9s \n', '#', 'z', 'LWC', 'R_eff');
fprintf(fileID, '%s %10s %9s %9s \n', '#', '(km)', '(g/m^3)', '(micron)');

% Lines run from the highest altitude to the surface. Precision: z to 0.1 mm
% and LWC to 1e-6 g m^-3. The old writer's 1 m / 1e-4 rounding put up to
% +/-0.4% scatter on tau and LWP for ~30 m synthetic slabs. libRadtran
% reads these columns in free format.
fprintf(fileID, '%12.7f %9.6f %9.4f \n', [z_lines'; lwc_lines'; re_lines']);
fclose(fileID);

end
