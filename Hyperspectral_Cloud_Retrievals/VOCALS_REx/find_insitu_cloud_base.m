%% Locate the measured cloud base of one in-situ vertical profile
%
% DESCRIPTION:
%   The profile finders select clouds using LWC and Nc thresholds that suit
%   the droplet retrievals, but those thresholds do not mark the physical cloud
%   base, and the adiabatic LWP depends on the cloud thickness squared. This
%   function uses the full flight time series to find the level where the cloud
%   water actually disappears. The profile selection itself is not changed; the
%   returned indices are used only for the adiabaticity estimate.
%
%   The base point is the first sample, moving away from the cloud, where the
%   droplet-mode LWC is below lwc_base_threshold_g_m3. At an adiabatic rate of
%   2 g/m^3 per km, 0.01 g/m^3 lies about 5 m above LWC = 0, which is the same
%   size as the height flown in one second, so that sample is a fair estimate
%   of where LWC reaches zero. The clear-air median of the droplet-mode LWC is
%   0 on every VOCALS-REx flight, so this threshold is above the noise floor.
%
% WHY THE RULE HAS THREE PARTS:
%
%   (1) Which end of the time series is the base depends on the flight
%       direction. On a climb the base comes first, so the search moves backward
%       in time. On a descent the base comes last, so it moves forward in time.
%
%   (2) DROPLET-MODE LWC, NOT TOTAL LWC. Drizzle falls out of the cloud base and
%       the 2DC keeps measuring it below the cloud. Searching on total LWC runs
%       as much as 170 m below the base on drizzling VOCALS-REx profiles, where
%       the droplet-mode LWC locates the base within about 15 m.
%
%   (3) THE PROFILE MAY ALREADY RUN PAST THE BASE. The profile finders end a
%       profile only when LWC AND Nc are BOTH below their thresholds, whereas a
%       profile starts only when both are above them. A descent can therefore
%       continue below the base while either quantity alone stays high (for
%       example drizzle keeping the total LWC above 0.05 g/m^3). When the
%       droplet-mode LWC at the base end of the profile is already below the
%       threshold, the function trims inward instead of walking outward. This
%       puts climbs and descents on the same footing.
%
% INPUTS:
%   altitude_m              - altitude for the WHOLE FLIGHT [m], 1 Hz
%   lwc_droplet_g_m3        - droplet-mode LWC for the WHOLE FLIGHT [g/m^3], 1 Hz.
%                             VOCALS-REx: lwc_CDP (no 2DC contribution)
%   i_start, i_end          - flight-series indices of the first and last sample
%                             of the profile, in time order
%   lwc_base_threshold_g_m3 - (optional, default 0.01) droplet-mode LWC below
%                             which the air counts as out of cloud [g/m^3]
%   n_consecutive           - (optional, default 1) how many consecutive samples
%                             must satisfy the condition before the search stops.
%                             Values above 1 protect against a one-second dropout
%                             in a broken cloud base being mistaken for the base.
%   max_walk_s              - (optional, default 180) longest outward search [s]
%
% OUTPUT:
%   base - structure with fields
%       .idx_profile            flight-series indices of the adiabaticity profile,
%                               in time order, ending at the base point
%       .idx_base               flight-series index of the base point
%       .is_ascending           true if the aircraft was climbing
%       .threshold_reached      false if no base was found; idx_profile is then the
%                               unmodified profile, so the base is not resolved
%       .n_walked_out           samples added below the original profile
%       .n_trimmed_in           samples removed from the original profile
%       .z_base_m               altitude of the base point [m]
%       .dz_below_original_m    lowest original sample minus lowest sample of the
%                               new profile [m]; negative when trimmed
%       .edge_gradient_g_m3_per_km  droplet-mode LWC gradient across the final
%                               step into the base point [g/m^3 per km]. Far above
%                               Gamma_ad (2 g/m^3 per km) this cannot be a parcel
%                               rising through the base, and more likely means the
%                               aircraft left through the side of the cloud. The
%                               true base is then below the base point found here.
%       .lwc_base_threshold_g_m3, .n_consecutive    settings used
%
% SEE ALSO:
%   compute_adiabaticity_Wood2005, find_verticalProfiles_VOCALS_REx_ver2
%
% By Andrew John Buggee
%%

function base = find_insitu_cloud_base(altitude_m, lwc_droplet_g_m3, i_start, i_end, ...
    lwc_base_threshold_g_m3, n_consecutive, max_walk_s)

if nargin < 5 || isempty(lwc_base_threshold_g_m3)
    lwc_base_threshold_g_m3 = 0.01;      % g/m^3
end

if nargin < 6 || isempty(n_consecutive)
    n_consecutive = 1;
end

if nargin < 7 || isempty(max_walk_s)
    max_walk_s = 180;                    % s
end

z_m = double(altitude_m(:))';            % m
lwc_g_m3 = double(lwc_droplet_g_m3(:))'; % g/m^3

% A profile needs a few samples left for its LWP integral to mean anything
min_samples_in_profile = 3;


% ------------------------------------------------------------------------
% Which end of the time series is the cloud base, and which way is "out"?
% ------------------------------------------------------------------------
is_ascending = mean(diff(z_m(i_start:i_end))) > 0;

if is_ascending
    i_edge = i_start;       % a climb meets cloud base first
    i_far_end = i_end;
    step_out = -1;          % leaving the cloud through its base means earlier times
else
    i_edge = i_end;         % a descent meets cloud base last
    i_far_end = i_start;
    step_out = +1;          % leaving the cloud through its base means later times
end


% ------------------------------------------------------------------------
% Search for the base point. Both branches leave i_base on the first sample
% below the threshold, so that LWC is about zero at the bottom of the profile.
% ------------------------------------------------------------------------
i_base = i_edge;
n_walked_out = 0;
n_trimmed_in = 0;
threshold_reached = false;
run_length = 0;

if lwc_g_m3(i_edge) >= lwc_base_threshold_g_m3

    % The profile ends inside cloud: walk outward, away from the cloud, until
    % the droplet-mode LWC has stayed below the threshold long enough.
    ii = i_edge;
    i_first_below = NaN;

    for kk = 1:max_walk_s

        ii = ii + step_out;

        if ii < 1 || ii > numel(lwc_g_m3)
            break       % ran off the end of the flight
        end

        if lwc_g_m3(ii) < lwc_base_threshold_g_m3

            if run_length == 0
                i_first_below = ii;
            end

            run_length = run_length + 1;

            if run_length >= n_consecutive
                threshold_reached = true;
                break
            end

        else
            run_length = 0;     % it was only a dropout, keep walking
        end

    end

    if threshold_reached
        i_base = i_first_below;
        n_walked_out = abs(i_base - i_edge);
    end

else

    % The profile already runs below the base: trim inward until the droplet-mode
    % LWC has stayed above the threshold long enough, keeping one sub-threshold
    % sample as the base point.
    ii = i_edge;
    i_first_above = NaN;

    while ii ~= i_far_end

        ii = ii - step_out;

        if lwc_g_m3(ii) >= lwc_base_threshold_g_m3

            if run_length == 0
                i_first_above = ii;
            end

            run_length = run_length + 1;

            if run_length >= n_consecutive
                threshold_reached = true;
                break
            end

        else
            run_length = 0;
        end

    end

    if threshold_reached
        i_base = i_first_above + step_out;
        n_trimmed_in = abs(i_base - i_edge);
    end

end


% ------------------------------------------------------------------------
% Build the profile that runs from the base point to the far end of the cloud.
% If trimming would leave almost nothing, the droplet mode never really formed
% a cloud here, so keep the original profile and report the base as not found.
% ------------------------------------------------------------------------
idx_profile = sort([i_base, i_far_end]);
idx_profile = idx_profile(1):idx_profile(2);

if numel(idx_profile) < min_samples_in_profile
    i_base = i_edge;
    idx_profile = i_start:i_end;
    n_walked_out = 0;
    n_trimmed_in = 0;
    threshold_reached = false;
end


% ------------------------------------------------------------------------
% Gradient across the final step into the base point: the sharp-edge diagnostic.
% The 0.5 m floor avoids dividing by zero if the aircraft momentarily levels off.
% ------------------------------------------------------------------------
i_inside = i_base - step_out;
dz_step_m = max(abs(z_m(i_inside) - z_m(i_base)), 0.5);                         % m
edge_gradient_g_m3_per_km = (lwc_g_m3(i_inside) - lwc_g_m3(i_base)) ./ dz_step_m .* 1e3;


base.idx_profile = idx_profile;
base.idx_base = i_base;
base.is_ascending = is_ascending;
base.threshold_reached = threshold_reached;
base.n_walked_out = n_walked_out;
base.n_trimmed_in = n_trimmed_in;
base.z_base_m = z_m(i_base);
base.dz_below_original_m = min(z_m(i_start:i_end)) - min(z_m(idx_profile));
base.edge_gradient_g_m3_per_km = edge_gradient_g_m3_per_km;
base.lwc_base_threshold_g_m3 = lwc_base_threshold_g_m3;
base.n_consecutive = n_consecutive;

end
