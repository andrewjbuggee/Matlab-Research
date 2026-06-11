function re_smooth_um = smooth_re_profile(re_um, z_m, window_m)
% Smooth an effective-radius profile with a vertical moving average.
%
% INPUTS:
%   re_um    - effective radius profile [microns], vector
%   z_m      - altitude of each level [meters], vector (same length as re_um)
%   window_m - full vertical width of the moving-average window [meters].
%              A value <= 0 (or []) returns the profile unchanged (the "raw"
%              case).
%
% OUTPUT:
%   re_smooth_um - smoothed effective radius profile [microns], with the same
%                  shape and altitude ordering as the input re_um.
%
% The moving average is computed in physical altitude units via movmean's
% 'SamplePoints' option, so the window represents the same vertical thickness
% regardless of how finely (or non-uniformly) a given profile is sampled.
% Window edges use the default 'shrink' behavior (the window shrinks toward
% cloud top/base rather than padding with zeros), so the boundary r_e values
% are not biased toward zero.
%
% This is the SAME smoothing used by both the preview script and the
% reflectance calculation, so what you see in the preview is exactly what is
% sent to libRadtran.
%
% By Andrew John Buggee

% --- preserve input orientation ---
was_row = isrow(re_um);
re_um   = re_um(:);
z_m     = z_m(:);

if numel(re_um) ~= numel(z_m)
    error([newline, 're_um and z_m must be the same length.', newline])
end

% --- raw profile requested ---
if isempty(window_m) || window_m <= 0
    re_smooth_um = re_um;
    if was_row, re_smooth_um = re_smooth_um'; end
    return
end

% --- sort ascending in altitude (SamplePoints expects increasing points) ---
[z_sorted, isort] = sort(z_m, 'ascend');
re_sorted = re_um(isort);

% --- break tied/duplicate altitudes so SamplePoints is strictly increasing
%     (a handful of in-situ profiles repeat an altitude level) ---
eps_z = 1e-3;   % meters
for ii = 2:numel(z_sorted)
    if z_sorted(ii) <= z_sorted(ii-1)
        z_sorted(ii) = z_sorted(ii-1) + eps_z;
    end
end

% --- moving average with the window defined in meters ---
re_sorted_smooth = movmean(re_sorted, window_m, 'SamplePoints', z_sorted, ...
    'Endpoints', 'shrink');

% --- restore the original ordering ---
re_smooth_um        = zeros(size(re_um));
re_smooth_um(isort) = re_sorted_smooth;

if was_row, re_smooth_um = re_smooth_um'; end

end
