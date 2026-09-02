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
%                  .tau            - cloud optical depth profile, computed from
%                                    Hansen & Travis (1974) Eq. (2.48) over the
%                                    measured size distribution at
%                                    tau_wavelength_nm. Stored in OPTICAL DEPTH
%                                    order: tau(1) = 0 at cloud top, tau(end) =
%                                    tau_c at cloud base. NOTE this is not the
%                                    same ordering as .altitude/.re/.lwc for
%                                    ascending profiles, which are stored in
%                                    sampling order -- flip those to pair them
%                                    with tau.

% By Andrew John Buggee

%%

function [vert_profs] = find_verticalProfiles_ORACLES(oracles, LWC_threshold, stop_at_max_lwc...
         ,Nc_threshold, which_computer)


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

% NOTE: significance_lvl and the nested chooseDistParam helper are no longer
% used. Both existed only to pick a fitted size-distribution shape parameter for
% the previous optical depth calculation, which has been replaced by the direct
% integral over the measured distribution. They are retained in case the fitted
% distributions are wanted for another purpose; nothing calls them today.
significance_lvl  = 0.1;    %#ok<NASGU> % 90% confidence for distribution fitting

% Wavelength at which the cloud optical depth is reported. Q_ext is evaluated
% here and the resulting tau is a monochromatic optical depth at this wavelength.
tau_wavelength_nm = 550;    % nm

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
    % ------------ Remove profiles that are not sampled over ocean ----------
    % -------------------------------------------------------------------------
    idx2delete = [];
    for nn = 1:length(vert_profs)

        % *** is the profile over ocean? ***
        isOcean = land_or_ocean(vert_profs(nn).latitude, vert_profs(nn).longitude, 20, false);
        if all(isOcean) == false
            idx2delete = [idx2delete, nn];
        end
    end
    if isempty(idx2delete)~= true
        vert_profs(idx2delete) = [];
    end


    

    % -------------------------------------------------------------------------
    % ------------- Remove near-duplicate profiles (multilayer edge cases) ----
    % -------------------------------------------------------------------------

    buffer_length = 4;
    idx2delete = [];

    for n1 = 1:length(vert_profs)
        for n2 = (n1+1):length(vert_profs)

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
    %
    % This is Hansen & Travis (1974) Equation (2.48) integrated over the depth of
    % the cloud, evaluated on the MEASURED droplet size distribution at every
    % height:
    %
    %     tau(z) = int_{z}^{z_top} [ int_r pi r^2 Q_ext(r,lambda) n(r,z') dr ] dz'
    %
    % The ORACLES probes report n(r,z) directly as Nd, the droplet concentration
    % per size bin (cm^-3), merged across the CAS, 2DS and HVPS-3 probes. Because
    % Nd is already a per-bin concentration rather than a density dN/dr, the
    % inner integral is a plain sum over bins with no bin-width factor:
    %
    %     beta(z) = pi * sum_i [ Q_ext(r_i) * r_i^2 * Nd_i(z) ]      [m^-1]
    %
    % NOTE ON THE PREVIOUS IMPLEMENTATION. This replaces a formulation that used
    %
    %     beta = pi * Q_ext_avg * r_eff^2 * total_Nc
    %
    % with Q_ext_avg from average_mie_over_size_distribution. That routine returns
    % the area-weighted efficiency of H&T (2.58), Q = k/G, which is defined to
    % pair with the geometric cross-section per unit volume, G = pi*N*<r^2> --
    % not with pi*N*r_eff^2. Since r_eff = <r^3>/<r^2> by H&T (2.53), Lyapunov's
    % moment inequality gives <r^2> <= r_eff^2 for every size distribution, so
    % that substitution biased tau high, by a median of about 35% across this
    % data set. See compare_tau_calculations_ORACLES.m for the full decomposition
    % and the numerical verification of each identity used here.
    %
    % STORAGE CONVENTION -- read this before using .tau alongside .re or .lwc.
    % tau is stored in OPTICAL DEPTH order: tau(1) = 0 at cloud top and
    % tau(end) = tau_c at cloud base, for every profile. The other profile fields
    % (altitude, re, lwc, Nd, ...) are stored in SAMPLING order, which runs
    % bottom-to-top for an ascending profile and top-to-bottom for a descending
    % one. So for ascending profiles tau is NOT element-aligned with the other
    % fields; consumers flip the others to match tau (see the fliplr calls in
    % ensemble_vertical_statistics_ORACLES.m). This preserves the convention the
    % previous implementation established and that downstream code relies on,
    % including every use of tau(end) as the total cloud optical depth.
    %
    % PERFORMANCE NOTE: Q_ext(r) depends only on radius and wavelength, not on
    % height or on the profile, so it is computed ONCE here for the whole data
    % set. The previous version ran a Mie calculation at every altitude of every
    % profile.

    if ~isempty(vert_profs)

        % Every profile shares the same probe size bins
        bin_centers_um = double(vert_profs(1).drop_radius_bin_center(:));    % um
        bin_edges_um   = double(vert_profs(1).drop_radius_bin_edges(:));     % um

        Q_ext_bin = binAveraged_Qext_forBins(bin_edges_um, bin_centers_um, ...
            tau_wavelength_nm, which_computer);

        r_bin_m2 = (bin_centers_um .* 1e-6).^2;         % um -> m, squared

    end

    for nn = 1:length(vert_profs)

        altitude_m = double(vert_profs(nn).altitude(:));      % m
        Nd_m3      = double(vert_profs(nn).Nd) .* 1e6;        % cm^-3 -> m^-3

        % --- H&T (2.48): extinction coefficient at every measured height ---
        beta = pi .* sum((Q_ext_bin .* r_bin_m2) .* Nd_m3, 1)';      % m^-1

        % --- Integrate downward from cloud top ---
        % Sorting by descending altitude makes the result independent of whether
        % the profile was sampled on a climb or a descent. The result is left in
        % this top-to-base order, which is the storage convention described above.
        [alt_desc, idx_desc] = sort(altitude_m, 'descend');

        depth_from_top = alt_desc(1) - alt_desc;              % m, 0 at cloud top
        tau_profile    = cumtrapz(depth_from_top, beta(idx_desc));

        if any(isnan(tau_profile))
            error([newline, 'NaN values found in optical depth for profile ', num2str(nn), newline])
        end

        vert_profs(nn).tau = reshape(tau_profile, 1, []);

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
            idx_start, idx_end)  %#ok<DEFNU>
        % Choose the best-fit size distribution model and return its scalar shape
        % parameter averaged over the cloud-top region [idx_start:idx_end].
        %
        % NO LONGER CALLED. The optical depth calculation that used this was
        % replaced by the direct H&T (2.48) integral over the measured size
        % distribution, which needs no fitted distribution. Kept for reference.

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
