%% Find ORACLES vertical cloud profiles (plane climbing or descending through cloud)

% Searches through ORACLES microphysics data to find contiguous segments
% where the P-3 aircraft flew vertically through a cloud layer. Computes
% liquid water path and cloud optical depth for each profile found.
%
% Analogous to find_verticalProfiles_VOCALS_REx_ver2.m, adapted for the
% ORACLES combined_microphysics netCDF dataset.

% INPUTS:
% -------
%   oracles        - structure returned by readORACLES.m
%
%   LWC_threshold  - (g/m^3) minimum LWC to define cloud boundaries.
%                    Uses the King hot-wire bulk_lwc, which is what the
%                    ORACLES data group used to define the cloud mask.
%                    The ORACLES dataset cloud definition is bulk_lwc >= 0.05.
%                    NOTE: do NOT use the distribution-computed lwc for
%                    this threshold — ~10% of in-cloud points have
%                    distribution lwc < 0.05 even when bulk_lwc >= 0.05,
%                    causing premature profile cutoffs and a low re bias.
%
%   stop_at_max_lwc - logical; if true, truncate each profile at peak LWC.
%
%   Nc_threshold   - (cm^-3) minimum total droplet number concentration.
%                    The ORACLES dataset cloud definition uses 10 cm^-3.
%
%   which_computer - string; used to locate Mie scattering tables.

% OUTPUTS:
% --------
%   vert_profs - array of structures, one per vertical profile found.
%                Each structure contains all oracles fields cropped to the
%                profile time window, plus:
%                  .LWC_threshold  - LWC threshold used [g/m^3]
%                  .Nc_threshold   - Nc threshold used [cm^-3]
%                  .lwp            - total liquid water path [g/m^2]
%                  .lwp_CAS        - CAS LWP (D <= 50 um) [g/m^2]
%                  .lwp_2DS_HVPS   - 2DS+HVPS LWP (D > 50 um) [g/m^2]
%                  .horz_dist      - cumulative horizontal distance [m]
%                  .cloud_depth    - geometric cloud depth [m]
%                  .slant_path     - slant path through cloud [m]
%                  .VR_zenith_angle - zenith angle of slant path [degrees]
%                  .tau            - cloud optical depth profile (0 at cloud top)

% By Andrew John Buggee

%%

function [vert_profs] = find_verticalProfiles_ORACLES(oracles, LWC_threshold, stop_at_max_lwc, Nc_threshold, which_computer)


% Grab all fieldnames for slicing
fields = fieldnames(oracles);


% -------------------------------------------------------------------------
% --------------------- Vertical profile requirements --------------------
% -------------------------------------------------------------------------
% Use bulk_lwc (King hot-wire) for all threshold comparisons.
% The ORACLES data group defined cloud as: N > 10 cm^-3 AND bulk_lwc > 0.05.
% After NaN->0 replacement in readORACLES, every non-zero bulk_lwc value
% is guaranteed to be >= 0.05, so it perfectly reproduces the cloud mask.
% The distribution-computed lwc has ~10% of in-cloud points below 0.05 due
% to variability in large-drop contributions, making it unsuitable here.
% -------------------------------------------------------------------------

% Compute instantaneous vertical velocity
dz_dt = diff(oracles.altitude) ./ diff(oracles.time);      % m/s

% Smooth with a sliding window to separate horizontal segments
n_window = 20;
dz_dt_mean = [movmean(dz_dt, n_window), 0];                % length N_time

% Vertical velocity threshold to distinguish ascending/descending
vertical_velocity_threshold = 2;    % m/s

% Index where plane is clearly ascending or descending
idx_dz_dt = abs(dz_dt_mean) > vertical_velocity_threshold;


% -------------------------------------------------------------------------
% ------------- Find transitions into and out of cloud --------------------
% -------------------------------------------------------------------------
% Use bulk_lwc for the cloud mask -- it never falls below 0.05 for in-cloud
% data, unlike the distribution-computed lwc.

idx_Nc_lwc_transition = [0, diff(oracles.total_Nc > Nc_threshold & ...
    oracles.bulk_lwc > LWC_threshold)];

% Starting indices of each potential cloud penetration
idx_1 = find(idx_Nc_lwc_transition == 1);

% Profile quality thresholds
length_threshold = 15;      % minimum number of data points in a profile
num_ba            = 20;      % points before/after profile to check for clear air
num_ba_long       = 30;      % extended clear-air check window

significance_lvl  = 0.1;    % 90% confidence for distribution fitting

profile_num = 0;


% -------------------------------------------------------------------------
% ---------------------- Step through cloud entries -----------------------
% -------------------------------------------------------------------------

for nn = 1:length(idx_1)

    % Find where both Nc and bulk_lwc next drop below their thresholds
    idx_cloud_boundary = find( ...
        oracles.total_Nc(idx_1(nn)+1:end)  < Nc_threshold & ...
        oracles.bulk_lwc(idx_1(nn)+1:end)  < LWC_threshold);

    if isempty(idx_cloud_boundary)
        continue
    end

    idx_cloud_boundary = idx_cloud_boundary(1) + idx_1(nn) - 1;

    % Was the plane ascending or descending throughout this cloud segment?
    ascend_or_descend_throughout_cloud = all(idx_dz_dt(idx_1(nn):idx_cloud_boundary));

    % Is the cloud segment long enough?
    meets_length_requirement = length(idx_1(nn):idx_cloud_boundary) > length_threshold;

    % Was the air before the profile below the thresholds?
    if (idx_1(nn) - num_ba) > 0
        before_profile_below_thresholds = ...
            median(oracles.total_Nc(idx_1(nn)-num_ba:idx_1(nn)-1)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(idx_1(nn)-num_ba:idx_1(nn)-1))    < LWC_threshold;
    else
        before_profile_below_thresholds = ...
            median(oracles.total_Nc(1:idx_1(nn)-1)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(1:idx_1(nn)-1))    < LWC_threshold;
    end

    % Was the air after the profile below the thresholds?
    if (idx_cloud_boundary + num_ba) < length(oracles.time)
        after_profile_below_thresholds = ...
            median(oracles.total_Nc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba))    < LWC_threshold;
    else
        after_profile_below_thresholds = ...
            median(oracles.total_Nc(idx_cloud_boundary+1:end)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(idx_cloud_boundary+1:end))    < LWC_threshold;
    end


    % -----------------------------------------------------------------
    % Case 1: clean single-layer profile
    % -----------------------------------------------------------------
    if ascend_or_descend_throughout_cloud && meets_length_requirement && ...
            before_profile_below_thresholds && after_profile_below_thresholds

        profile_num = profile_num + 1;
        vert_profs(profile_num) = oracles;
        vert_profs(profile_num) = cropProfileFields(vert_profs(profile_num), fields, ...
            idx_1(nn), idx_cloud_boundary, idx_dz_dt);


        % -----------------------------------------------------------------
        % Case 2: possible multi-layer cloud -- gap after profile
        % -----------------------------------------------------------------
    elseif ascend_or_descend_throughout_cloud && meets_length_requirement && ...
            before_profile_below_thresholds && ~after_profile_below_thresholds

        % Find where the next cloud layer begins
        idx_cloud_boundary2 = find( ...
            oracles.total_Nc(idx_cloud_boundary+1:end) > Nc_threshold & ...
            oracles.bulk_lwc(idx_cloud_boundary+1:end) > LWC_threshold);

        if isempty(idx_cloud_boundary2)
            continue
        end
        idx_cloud_boundary2 = idx_cloud_boundary2(1) + idx_cloud_boundary;

        length_between_multi_layers = idx_cloud_boundary2 - idx_cloud_boundary - 1;

        % Find where the second layer ends
        idx_cloud_boundary3 = find( ...
            oracles.total_Nc(idx_cloud_boundary2+1:end) < Nc_threshold & ...
            oracles.bulk_lwc(idx_cloud_boundary2+1:end) < LWC_threshold);

        if isempty(idx_cloud_boundary3)
            continue
        end
        idx_cloud_boundary3 = idx_cloud_boundary3(1) + idx_cloud_boundary2 - 1;

        if (idx_cloud_boundary + num_ba_long) >= length(oracles.time)
            continue
        end
        after_profile_below_thresholds_long = ...
            median(oracles.total_Nc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba_long)) < Nc_threshold & ...
            median(oracles.bulk_lwc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba_long)) < LWC_threshold;

        if (idx_cloud_boundary3 + num_ba) >= length(oracles.time)
            continue
        end
        after_profile_below_thresholds2 = ...
            median(oracles.total_Nc(idx_cloud_boundary3:idx_cloud_boundary3+num_ba)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(idx_cloud_boundary3:idx_cloud_boundary3+num_ba))    < LWC_threshold;

        if length_between_multi_layers < 5 && after_profile_below_thresholds_long && ...
                after_profile_below_thresholds2

            profile_num = profile_num + 1;
            vert_profs(profile_num) = oracles;
            vert_profs(profile_num) = cropProfileFields(vert_profs(profile_num), fields, ...
                idx_1(nn), idx_cloud_boundary3, idx_dz_dt);
        end


        % -----------------------------------------------------------------
        % Case 3: possible multi-layer cloud -- gap before profile
        % -----------------------------------------------------------------
    elseif ascend_or_descend_throughout_cloud && meets_length_requirement && ...
            ~before_profile_below_thresholds && after_profile_below_thresholds

        % Find where the preceding cloud layer ended
        idx_cloud_boundary0 = find( ...
            oracles.total_Nc(1:idx_1(nn)-1) > Nc_threshold & ...
            oracles.bulk_lwc(1:idx_1(nn)-1) > LWC_threshold);

        if isempty(idx_cloud_boundary0)
            continue
        end
        idx_cloud_boundary0 = idx_cloud_boundary0(end);

        length_between_multi_layers = idx_1(nn) - idx_cloud_boundary0;

        % Find where the preceding layer began
        idx_cloud_boundary_minus1 = find( ...
            oracles.total_Nc(1:idx_cloud_boundary0-1) < Nc_threshold & ...
            oracles.bulk_lwc(1:idx_cloud_boundary0-1) < LWC_threshold);

        if isempty(idx_cloud_boundary_minus1)
            continue
        end
        idx_cloud_boundary_minus1 = idx_cloud_boundary_minus1(end);

        length_of_second_layer = idx_cloud_boundary0 - idx_cloud_boundary_minus1;

        if (idx_1(nn) - num_ba_long) < 1
            continue
        end
        before_profile_below_thresholds_long = ...
            median(oracles.total_Nc(idx_1(nn)-num_ba_long:idx_1(nn)-1)) < Nc_threshold & ...
            median(oracles.bulk_lwc(idx_1(nn)-num_ba_long:idx_1(nn)-1)) < LWC_threshold;

        if (idx_cloud_boundary_minus1 - num_ba_long) < 1
            continue
        end
        before_profile_below_thresholds2 = ...
            median(oracles.total_Nc(idx_cloud_boundary_minus1-num_ba_long:idx_cloud_boundary_minus1)) < Nc_threshold & ...
            mean(oracles.bulk_lwc(idx_cloud_boundary_minus1-num_ba_long:idx_cloud_boundary_minus1))    < LWC_threshold;

        if length_between_multi_layers < 5 && before_profile_below_thresholds_long && ...
                before_profile_below_thresholds2

            profile_num = profile_num + 1;
            vert_profs(profile_num) = oracles;
            vert_profs(profile_num) = cropProfileFields(vert_profs(profile_num), fields, ...
                idx_cloud_boundary_minus1, idx_cloud_boundary, idx_dz_dt);

        elseif length_between_multi_layers > length_of_second_layer && ...
                before_profile_below_thresholds_long && before_profile_below_thresholds2

            profile_num = profile_num + 1;
            vert_profs(profile_num) = oracles;
            vert_profs(profile_num) = cropProfileFields(vert_profs(profile_num), fields, ...
                idx_1(nn), idx_cloud_boundary, idx_dz_dt);
        end

    end


end

% If there are no vertical profiles found, terminate the function
if exist("vert_profs", "var") == 0

    % end the function with a statement that no vert profiles were found
    vert_profs = [];

else



    % -------------------------------------------------------------------------
    % ------------- Remove near-duplicate profiles (multilayer edge cases) ----
    % -------------------------------------------------------------------------

    buffer_length = 4;
    idx2delete = [];

    for n1 = 1:profile_num
        for n2 = (n1+1):profile_num

            if (vert_profs(n2).time(1) - vert_profs(n1).time(1)) <= buffer_length || ...
                    (vert_profs(n2).time(end) - vert_profs(n1).time(end)) <= buffer_length

                if length(vert_profs(n2).time) > length(vert_profs(n1).time)
                    idx2delete = [idx2delete, n1];
                elseif length(vert_profs(n2).time) < length(vert_profs(n1).time)
                    idx2delete = [idx2delete, n2];
                else
                    idx2delete = [idx2delete, n2];
                end

            end

        end
    end

    if isempty(idx2delete)~= true
        vert_profs(idx2delete) = [];
    end


    % -------------------------------------------------------------------------
    % --------- Remove profiles with too many below-threshold Nc points -------
    % -------------------------------------------------------------------------

    idx2delete = [];
    for nn = 1:length(vert_profs)
        if sum(vert_profs(nn).total_Nc < Nc_threshold) > floor(length(vert_profs(nn).total_Nc)/2)
            idx2delete = [idx2delete, nn];
        end
    end
    if isempty(idx2delete)~= true
        vert_profs(idx2delete) = [];
    end


    % -------------------------------------------------------------------------
    % -------------- Optionally truncate at peak LWC --------------------------
    % -------------------------------------------------------------------------

    if stop_at_max_lwc
        error([newline, 'The code to cut profiles at max LWC has not been tested.', newline])
    end


    % -------------------------------------------------------------------------
    % ------------- Store threshold values in each profile --------------------
    % -------------------------------------------------------------------------

    for nn = 1:length(vert_profs)
        vert_profs(nn).LWC_threshold = LWC_threshold;   % g/m^3
        vert_profs(nn).Nc_threshold  = Nc_threshold;     % cm^-3
    end


    % -------------------------------------------------------------------------
    % ---------------------- Compute Liquid Water Path ------------------------
    % -------------------------------------------------------------------------

    for nn = 1:length(vert_profs)

        dz_dt_prof = diff(vert_profs(nn).altitude) ./ diff(vert_profs(nn).time);

        if mean(dz_dt_prof) > 0
            sign_factor = 1;        % ascending
        else
            sign_factor = -1;       % descending
        end

        vert_profs(nn).lwp          = sign_factor * trapz(vert_profs(nn).altitude, vert_profs(nn).lwc);       % g/m^2
        vert_profs(nn).lwp_CAS      = sign_factor * trapz(vert_profs(nn).altitude, vert_profs(nn).cwc);       % g/m^2
        vert_profs(nn).lwp_2DS_HVPS = sign_factor * trapz(vert_profs(nn).altitude, vert_profs(nn).rwc);       % g/m^2

    end


    % -------------------------------------------------------------------------
    % -------------- Compute horizontal distance and slant path ---------------
    % -------------------------------------------------------------------------

    wgs84 = wgs84Ellipsoid("m");

    for nn = 1:length(vert_profs)

        horz_distance_travelled = zeros(1, length(vert_profs(nn).latitude));

        for xx = 2:length(vert_profs(nn).latitude)
            horz_distance_travelled(xx) = distance( ...
                vert_profs(nn).latitude(1),  vert_profs(nn).longitude(1), ...
                vert_profs(nn).latitude(xx), vert_profs(nn).longitude(xx), wgs84);
        end

        vert_profs(nn).horz_dist       = horz_distance_travelled;
        vert_profs(nn).cloud_depth     = max(vert_profs(nn).altitude) - min(vert_profs(nn).altitude);
        vert_profs(nn).slant_path      = sqrt(vert_profs(nn).horz_dist(end)^2 + vert_profs(nn).cloud_depth^2);
        vert_profs(nn).VR_zenith_angle = atand(vert_profs(nn).horz_dist(end) / vert_profs(nn).cloud_depth);

    end


    % -------------------------------------------------------------------------
    % ---------------------- Compute optical depth ----------------------------
    % -------------------------------------------------------------------------
    % Optical depth is defined as 0 at cloud top, increasing toward cloud bottom.
    % Uses the ORACLES pre-computed effective radius (re), which incorporates
    % all three probes (CAS, 2DS, HVPS-3).
    %
    % PERFORMANCE NOTE: find_bestFitDist_dropDist is called ONCE per profile,
    % outside the inner altitude loop. The distribution fit over the full
    % profile is constant, so recomputing it at every altitude step is wasteful.
    % Only average_mie_over_size_distribution (which depends on the cumulative
    % re profile) is called inside the inner loop.

    for nn = 1:length(vert_profs)

        vector_length = length(vert_profs(nn).altitude);
        vert_profs(nn).tau = zeros(1, vector_length - 1);

        dz_dt_prof   = diff(vert_profs(nn).altitude) ./ diff(vert_profs(nn).time);
        is_ascending = mean(dz_dt_prof) > 0;

        % ------------------------------------------------------------------
        % Fit size distribution ONCE per profile (outside the inner loop).
        % The result is the same for every altitude step in this profile.
        % ------------------------------------------------------------------
        [normFit, logNormFit, gammaFit] = find_bestFitDist_dropDist( ...
            vert_profs(nn).Nd, vert_profs(nn).drop_radius_bin_edges, ...
            vert_profs(nn).drop_radius_bin_center, significance_lvl);

        % Average shape parameter over the cloud-top region, which dominates
        % the optical depth integral regardless of ascent/descent direction
        if is_ascending
            idx_shape_start = round(2*vector_length/3);
            idx_shape_end   = vector_length;
        else
            idx_shape_start = 1;
            idx_shape_end   = round(vector_length/3);
        end

        % Select distribution model and pre-compute the scalar shape parameter
        [dist_2model, shape_param] = chooseDistParam( ...
            normFit, logNormFit, gammaFit, idx_shape_start, idx_shape_end);


        if is_ascending
            % --- Ascending: cloud bottom sampled first, cloud top last ---

            for ii = 1:vector_length - 1
                
                % for debugging
                %disp([newline, 'iteration: ', num2str(ii), newline])

                % Cumulative integration from cloud top downward to level ii
                re_meters       = vert_profs(nn).re(vector_length-ii:vector_length) ./ 1e6;   % m
                total_Nc_meters = vert_profs(nn).total_Nc(vector_length-ii:vector_length) .* 1e6; % m^-3
                alt_from_top    = vert_profs(nn).altitude(end) - ...
                    vert_profs(nn).altitude(vector_length-ii:vector_length);   % m

                re_meters(re_meters == 0) = 0.009e-6;   % guard against zero re at edges

                % Build constant-valued distribution vector at the correct length
                distribution_dist = repmat(shape_param, 1, length(re_meters));
                  
                [~, Qe_avg, ~] = average_mie_over_size_distribution( ...
                    re_meters .* 1e6, distribution_dist, 550, 'water', dist_2model, which_computer, ii);

                vert_profs(nn).tau(ii) = pi * trapz(fliplr(alt_from_top), ...
                    fliplr(Qe_avg .* re_meters.^2 .* total_Nc_meters));

            end

        else
            % --- Descending: cloud top sampled first, cloud bottom last ---

            for ii = 1:vector_length - 1

                re_meters       = vert_profs(nn).re(1:ii+1) ./ 1e6;   % m
                total_Nc_meters = vert_profs(nn).total_Nc(1:ii+1) .* 1e6; % m^-3
                alt_from_top    = vert_profs(nn).altitude(1) - vert_profs(nn).altitude(1:ii+1);   % m

                re_meters(re_meters == 0) = 0.009e-6;

                distribution_dist = repmat(shape_param, 1, length(re_meters));

                [~, Qe_avg, ~] = average_mie_over_size_distribution( ...
                    re_meters .* 1e6, distribution_dist, 550, 'water', dist_2model, which_computer, ii);

                vert_profs(nn).tau(ii) = pi * trapz(alt_from_top, ...
                    Qe_avg(:,end) .* re_meters.^2 .* total_Nc_meters);

            end

        end

        if sum(isnan(vert_profs(nn).tau)) > 0
            error([newline, 'NaN values found in optical depth for profile ', num2str(nn), newline])
        end

        % Prepend zero: tau = 0 at cloud top
        vert_profs(nn).tau = [0, vert_profs(nn).tau];

    end


end


% =========================================================================
% =================== Local helper functions ==============================
% =========================================================================

    function prof = cropProfileFields(prof, fields, idx_start, idx_end, idx_dz_dt)
        % Crop all time-varying fields of prof to the range [idx_start, idx_end]

        N_time = length(idx_dz_dt);

        for ff = 1:length(fields)

            n_elem = numel(prof.(fields{ff}));

            if n_elem == N_time
                % 1D time-varying field: crop to row vector
                prof.(fields{ff}) = reshape(prof.(fields{ff})(idx_start:idx_end), 1, []);

            elseif n_elem > N_time
                % 2D field (bins x time): crop along time dimension
                prof.(fields{ff}) = prof.(fields{ff})(:, idx_start:idx_end);

                % else: non-time field (bin edges, scalars, datetime) -- keep as-is
            end

        end

    end


    function [dist_2model, shape_param] = chooseDistParam(normFit, logNormFit, gammaFit, ...
            idx_start, idx_end)
        % Choose the best-fit size distribution model and return its scalar shape
        % parameter averaged over the cloud-top region [idx_start:idx_end].

        h_norm    = sum(normFit.h_test,    "omitnan");
        h_lognorm = sum(logNormFit.h_test, "omitnan");
        h_gamma   = sum(gammaFit.h_test,   "omitnan");

        if h_lognorm < h_norm && h_lognorm < h_gamma
            dist_2model = 'lognormal';
            shape_param = mean(logNormFit.std(idx_start:idx_end), "omitnan");

        elseif h_gamma <= h_norm && h_gamma <= h_lognorm
            dist_2model = 'gamma';
            shape_param = mean(gammaFit.alpha(idx_start:idx_end), "omitnan");

        else
            % Default / tie: gamma with alpha = 10
            dist_2model = 'gamma';
            if isfield(gammaFit, 'alpha') && ~isempty(gammaFit.alpha)
                shape_param = mean(gammaFit.alpha(idx_start:idx_end), "omitnan");
            else
                shape_param = 10;
            end

        end

        % Guard against NaN or non-positive shape parameter (thin cloud edge case)
        if isnan(shape_param) || shape_param <= 0
            shape_param = 10;
        end

    end



end
