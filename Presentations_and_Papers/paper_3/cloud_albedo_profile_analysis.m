%% Cloud Albedo vs. Vertical Droplet Profile Analysis
%
% Investigates how the vertical profile of droplet effective radius impacts
% broadband cloud albedo, motivated by Borg & Bennartz (2007, GRL).
%
% Two main comparisons are performed:
%
%   Part 1 — Fixed LWP (following Borg & Bennartz 2007)
%     For each VOCALS-REx in-situ profile, three idealized clouds are
%     constructed with identical LWP and cloud boundaries but different r_e
%     profile shapes:
%       (a) True in-situ measured r_e profile
%       (b) Smooth monotonic adiabatic profile:  r_e ∝ (z – z_base)^(1/3)
%       (c) Constant r_e profile (LWC-weighted mean)
%
%   Part 2 — Retrieval comparison (same cloud boundaries and tau)
%     For each profile, three clouds are constructed representing what
%     different retrievals would "see":
%       (a) True in-situ profile
%       (b) Linear profile retrieval (r_top, r_bot) — as in Buggee & Pilewskie (2026)
%       (c) Two-band LUT (TBLUT) retrieval — constant r_e, reflectance-weighted
%           toward cloud top (single r_eff, same tau_c)
%
%   Part 3 — Cloud-top droplet gradient analysis
%     Find profiles with the steepest r_e change near cloud top.
%     Connect those gradients to:
%       * Cloud-top entrainment instability (CTEI, κ parameter)
%       * Albedo bias between retrieval cases
%
% Broadband (250–4000 nm) solar albedo is computed via libRadtran (DISORT,
% 16 streams), using the same atmospheric settings as the HySICS training
% simulations in hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver3.m.
%
% VOCALS-REx and ORACLES ensemble profiles are supported via the 'campaign'
% switch in the data-loading section.
%
% By Andrew John Buggee

clear variables

%% -------------------------------------------------------------------------
%  Paths and setup
%  -------------------------------------------------------------------------

addLibRadTran_paths;

which_computer = whatComputer();

% ---- libRadtran folder paths ----
if strcmp(which_computer, 'andrewbuggee')

    folder_paths = define_folderPaths_for_HySICS(5);

    foldername_vocals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_oracles = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/';

elseif strcmp(which_computer, 'anbu8374')

    folder_paths.which_computer           = which_computer;
    folder_paths.libRadtran_inp           = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/';
    folder_paths.libRadtran_data          = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';
    folder_paths.libRadtran_water_cloud_files = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';
    folder_paths.atm_folder_path          = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';

    foldername_vocals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_oracles = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/';

elseif strcmp(which_computer, 'curc')

    folder_paths.which_computer           = which_computer;
    folder_paths.libRadtran_inp           = '/scratch/alpine/anbu8374/HySICS/INP_OUT/';
    folder_paths.libRadtran_data          = '/projects/anbu8374/software/libRadtran-2.0.5/data/';
    folder_paths.libRadtran_water_cloud_files = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';
    folder_paths.atm_folder_path          = '/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod/';

    foldername_vocals = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
        'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_oracles = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/', ...
        'ORACLES/oracles_data/'];

    foldername_save = '/projects/anbu8374/Matlab-Research/Presentations_and_Papers/paper_3/';

end


%% -------------------------------------------------------------------------
%  User settings
%  -------------------------------------------------------------------------

% Choose campaign: 'vocalsRex' or 'oracles'
campaign_name = 'vocalsRex';

% Solar zenith angle [degrees] — representative value for VOCALS-REx region
sza = 30;

% Solar azimuth angle [degrees] — 0 = due South (libRadtran convention)
phi0 = 0;

% Day of year (affects Earth–Sun distance correction)
% 316 = 12 November (VOCALS-REx period); change as needed
day_of_year = 316;

% Surface albedo (ocean)
surface_albedo = 0.04;

% Drizzle LWP threshold used to flag precipitating profiles [g/m²]
drizzle_LWP_threshold = 5;

% Fraction of cloud depth from cloud top used to compute cloud-top r_e gradient
cloud_top_fraction = 0.15;   % top 15% of cloud depth


%% -------------------------------------------------------------------------
%  Physical constants
%  -------------------------------------------------------------------------

g     = 9.81;          % m/s²
c_pd  = 1005;          % J/(kg·K)
L_v   = 2.501e6;       % J/kg
R_d   = 287.04;        % J/(kg·K)
R_v   = 461.5;         % J/(kg·K)
epsilon = R_d / R_v;   % ≈ 0.622
rho_L = 1000;          % kg/m³  (liquid water density)


%% -------------------------------------------------------------------------
%  Load ensemble profiles
%  -------------------------------------------------------------------------

if strcmp(campaign_name, 'vocalsRex')

    % Most recent VOCALS-REx non-precipitating ensemble file
    vocals_filename = ['ensemble_profiles_without_precip_from_14_files_LWC-threshold',...
        '_0.03_Nc-threshold_25_drizzleLWP-threshold_5_04-Dec-2025.mat'];

    if ~exist([foldername_vocals, vocals_filename], 'file')
        % Fall back: find any ensemble_profiles file
        d = dir([foldername_vocals, 'ensemble_profiles_without_precip*.mat']);
        if isempty(d)
            error('No ensemble profiles .mat file found in:\n  %s', foldername_vocals)
        end
        vocals_filename = d(end).name;
    end

    disp(['Loading VOCALS-REx: ', vocals_filename])
    tmp = load([foldername_vocals, vocals_filename], 'ensemble_profiles');
    ensemble_profiles = tmp.ensemble_profiles;

elseif strcmp(campaign_name, 'oracles')

    d = dir([foldername_oracles, 'ensemble_profiles*.mat']);
    if isempty(d)
        error('No ensemble profiles .mat file found in:\n  %s', foldername_oracles)
    end
    oracles_filename = d(end).name;
    disp(['Loading ORACLES: ', oracles_filename])
    tmp = load([foldername_oracles, oracles_filename], 'ensemble_profiles');
    ensemble_profiles = tmp.ensemble_profiles;

end

N_profiles = length(ensemble_profiles);
fprintf('Loaded %d vertical profiles from %s.\n', N_profiles, campaign_name)


%% -------------------------------------------------------------------------
%  Pre-screen profiles
%  Exclude precipitating profiles and any profile with too few levels
%  -------------------------------------------------------------------------

is_drizzle  = false(1, N_profiles);
is_valid    = true(1, N_profiles);

for nn = 1:N_profiles
    prof = ensemble_profiles{nn};

    % Drizzle flag
    if isfield(prof, 'lwp_2DS_HVPS')
        is_drizzle(nn) = prof.lwp_2DS_HVPS >= drizzle_LWP_threshold;

    elseif isfield(prof, 'lwp_2DC')
        is_drizzle(nn) = prof.lwp_2DC >= drizzle_LWP_threshold;
    end

    % Orient the profile and check minimum length
    [re_p, lwc_p, ~, alt_p, ~, ~, ~] = orientProfile(prof);

    if length(re_p) < 5 || any(isnan(re_p)) || any(re_p <= 0)
        is_valid(nn) = false;
    end

    if any(lwc_p < 0)
        is_valid(nn) = false;
    end
end

idx_use = find(is_valid & ~is_drizzle);
fprintf('Profiles used (non-drizzling, valid): %d of %d\n', length(idx_use), N_profiles)


%% =========================================================================
%  PART 1 — Fixed-LWP albedo comparison (Borg & Bernard 2007 approach)
%  =========================================================================
%
%  For each profile the same measured LWC distribution (→ same LWP) is
%  kept.  Three effective-radius profile shapes are imposed:
%    (a) in-situ measured r_e
%    (b) adiabatic shape  r_e ∝ (z – z_base)^(1/3), scaled to same LWP-tau
%    (c) constant r_e    = LWC-weighted mean of in-situ r_e
%  =========================================================================

fprintf('\n--- Part 1: Fixed-LWP albedo comparison ---\n')

N_use = length(idx_use);

% Pre-allocate albedo storage
albedo_fixedLWP.insitu    = NaN(1, N_use);
albedo_fixedLWP.adiabatic = NaN(1, N_use);
albedo_fixedLWP.constant  = NaN(1, N_use);

% Profile-level diagnostics
LWP_all           = NaN(1, N_use);
tau_insitu        = NaN(1, N_use);
tau_adiabatic     = NaN(1, N_use);
tau_constant      = NaN(1, N_use);
cloudTop_re_grad  = NaN(1, N_use);   % d(r_e)/d(z) near cloud top [µm/m]
re_top_insitu     = NaN(1, N_use);   % r_e at cloud top [µm]
re_base_insitu    = NaN(1, N_use);   % r_e at cloud base [µm]
cloud_depth_all   = NaN(1, N_use);   % geometric cloud depth [m]


for kk = 1:N_use

    nn = idx_use(kk);
    prof = ensemble_profiles{nn};

    % ---- Orient to cloud-base-first ----
    [re_p, lwc_p, Nc_p, alt_p, T_K_p, P_Pa_p, ~] = orientProfile(prof);
    % alt_p(1) = cloud base, alt_p(end) = cloud top

    % Remove any remaining zero/negative droplet levels
    valid_lev = (re_p > 0.5) & (lwc_p > 0) & ~isnan(re_p);
    if sum(valid_lev) < 5
        continue
    end
    re_p   = re_p(valid_lev);
    lwc_p  = lwc_p(valid_lev);
    Nc_p   = Nc_p(valid_lev);
    alt_p  = alt_p(valid_lev);
    T_K_p  = T_K_p(valid_lev);
    P_Pa_p = P_Pa_p(valid_lev);

    % Convert altitude to km for wc-file writer
    z_km   = alt_p(:) ./ 1e3;         % km MSL, column vector
    re_col = re_p(:);                  % µm
    lwc_col = lwc_p(:);                % g/m³

    % ---- Scalar cloud properties ----
    H_total = abs(alt_p(end) - alt_p(1));        % geometric depth [m]
    LWP     = trapz(alt_p, lwc_p);               % g/m²
    LWP_all(kk) = LWP;
    cloud_depth_all(kk) = H_total;

    % ---- Cloud-top r_e gradient ----
    % Use top fraction of cloud depth
    z_from_top  = alt_p(end) - alt_p;
    top_idx     = z_from_top <= cloud_top_fraction * H_total;
    if sum(top_idx) >= 2
        dre_dz = (re_p(end) - re_p(find(top_idx, 1, 'first'))) / ...
            (alt_p(end) - alt_p(find(top_idx, 1, 'first')));
        cloudTop_re_grad(kk) = dre_dz;   % µm/m (negative = decreasing toward top)
    end

    re_top_insitu(kk)  = re_p(end);   % cloud top value
    re_base_insitu(kk) = re_p(1);     % cloud base value

    % ---- Cloud-top optical depth ----
    % tau_c ≈ (3/2) * integral(LWC/r_e) dz / rho_L
    %       [with LWC in g/m³, r_e in µm, rho_L=1e6 g/m³, dz in m]
    %  tau = (3/2) * integral(LWC * 1e-3 / (r_e * 1e-6)) dz
    %       = (3/2) * integral(LWC / r_e) * 1e3 dz  [m⁻¹ integrated over m]
    % tau_insitu(kk) = (3/2) * trapz(alt_p, (lwc_p ./ re_p) * 1e3 / rho_L);
    tau_insitu(kk) = max(prof.tau);


    % =======================================================================
    % (b)  Adiabatic r_e profile with same measured LWC distribution
    % =======================================================================
    % Build an adiabatic-shaped r_e profile that preserves the measured LWC:
    %   r_e_ad(z) = A * (z - z_base)^(1/3)
    %   The constant A is set so that the mean (LWC-weighted) r_e matches
    %   the in-situ LWC-weighted mean, preserving the total optical depth.

    z_from_base = max(alt_p - alt_p(1), 0);   % [m]

    % Adiabatic shape  (avoid z=base singularity by using a small offset)
    re_shape_ad = (z_from_base + 1).^(1/3);   % [m^(1/3)], offset of 1 m to avoid 0

    % Scale to match LWC-weighted mean r_e of in-situ profile
    lwc_weights     = lwc_p / sum(lwc_p);
    re_mean_insitu  = sum(lwc_weights .* re_p);
    re_mean_shape   = sum(lwc_weights .* re_shape_ad);
    scale_ad        = re_mean_insitu / re_mean_shape;
    re_adiabatic    = re_shape_ad * scale_ad;   % µm — same LWC-weighted mean

    tau_adiabatic(kk) = (3/2) * trapz(alt_p, (lwc_p ./ re_adiabatic) * 1e3 / rho_L);


    % =======================================================================
    % (c)  Constant r_e profile with same measured LWC
    % =======================================================================
    % r_e_const = LWC-weighted mean of in-situ r_e → same LWC-weighted mean
    % This is equivalent to what a column-mean retrieval would assume.

    re_constant = repmat(re_mean_insitu, size(re_p));

    tau_constant(kk) = (3/2) * trapz(alt_p, (lwc_p ./ re_constant) * 1e3 / rho_L);


    % =======================================================================
    % Run libRadtran for all three profile types
    % =======================================================================

    % Date and time metadata for wc-file naming
    if isfield(prof, 'dateOfFlight')
        date_str = prof.dateOfFlight;
    else
        date_str = '2008-11-09';
    end

    if isfield(prof, 'time_utc')
        time_val = mean(prof.time_utc(~isnan(prof.time_utc)), 'omitnan');
    else
        time_val = 120000;
    end

    % Ensure z increases from cloud top to bottom for wc file
    % write_wc_file_from_in_situ expects top-to-bottom ordering
    if z_km(2) > z_km(1)
        % currently bottom-to-top; flip
        z_wc = flipud(z_km);
        lwc_wc_in    = flipud(lwc_col);
        re_wc_insitu = flipud(re_col);
        re_wc_ad     = flipud(re_adiabatic(:));
        re_wc_const  = flipud(re_constant(:));
    else
        z_wc          = z_km;
        lwc_wc_in     = lwc_col;
        re_wc_insitu  = re_col;
        re_wc_ad      = re_adiabatic(:);
        re_wc_const   = re_constant(:);
    end

    % Build libRadtran settings structure for this profile
    RT.sza            = sza;
    RT.phi0           = phi0;
    RT.day_of_year    = day_of_year;
    RT.surface_albedo = surface_albedo;
    RT.atm_file       = 'afglus.dat';
    RT.source_file    = 'binned_fs_hybrid_reference_spectrum_c2022-11-30_with_unc.dat';     % Need full extension reference spectrum

    % -- (a) In-situ --
    try
        albedo_fixedLWP.insitu(kk) = compute_broadband_albedo(...
            re_wc_insitu, lwc_wc_in, z_wc, ...
            RT, folder_paths, kk, 'insitu', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (in-situ, Part1): %s', kk, ME.message)
    end

    % -- (b) Adiabatic --
    try
        albedo_fixedLWP.adiabatic(kk) = compute_broadband_albedo(...
            re_wc_ad, lwc_wc_in, z_wc, ...
            RT, folder_paths, kk, 'adiabatic', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (adiabatic, Part1): %s', kk, ME.message)
    end

    % -- (c) Constant --
    try
        albedo_fixedLWP.constant(kk) = compute_broadband_albedo(...
            re_wc_const, lwc_wc_in, z_wc, ...
            RT, folder_paths, kk, 'constant', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (constant, Part1): %s', kk, ME.message)
    end

    fprintf('  Part 1: profile %d/%d  |  albedo: insitu=%.4f  adiab=%.4f  const=%.4f\n', ...
        kk, N_use, albedo_fixedLWP.insitu(kk), ...
        albedo_fixedLWP.adiabatic(kk), albedo_fixedLWP.constant(kk))

end % kk loop Part 1


%% =========================================================================
%  PART 2 — Retrieval comparison (same tau_c and cloud boundaries)
%  =========================================================================
%
%  Three "retrieved" clouds are constructed from each in-situ profile,
%  all with the same optical depth and boundaries but different r_e profiles:
%
%    (a) True in-situ profile  — ground truth
%    (b) Linear profile retrieval: r_e varies linearly from r_bot (base)
%        to r_top (top), parameterised as in Buggee & Pilewskie (2026).
%        r_top = r_e at cloud top (in-situ), r_bot = r_e at cloud base.
%        LWC is set from tau_c via LWC = r_e * tau_step * (2/3) * rho_L.
%    (c) TBLUT-style: constant r_e = reflectance-weighted effective radius
%        (weighted toward cloud top, mimicking 2-band retrieval sensitivity)
%        Same tau_c as in-situ.
%  =========================================================================

fprintf('\n--- Part 2: Retrieval comparison (same tau) ---\n')

albedo_retrieval.insitu   = NaN(1, N_use);
albedo_retrieval.profile  = NaN(1, N_use);
albedo_retrieval.tblut    = NaN(1, N_use);

re_tblut_all     = NaN(1, N_use);   % store the TBLUT r_eff for each profile
re_linear_top    = NaN(1, N_use);   % r_top from linear profile retrieval
re_linear_bot    = NaN(1, N_use);   % r_bot from linear profile retrieval


for kk = 1:N_use

    nn  = idx_use(kk);
    prof = ensemble_profiles{nn};

    [re_p, lwc_p, Nc_p, alt_p, T_K_p, P_Pa_p, ~] = orientProfile(prof);

    valid_lev = (re_p > 0.5) & (lwc_p > 0) & ~isnan(re_p);
    if sum(valid_lev) < 5
        continue
    end
    re_p   = re_p(valid_lev);
    lwc_p  = lwc_p(valid_lev);
    alt_p  = alt_p(valid_lev);

    z_km   = alt_p(:) ./ 1e3;   % km, bottom-to-top
    N_lev  = length(alt_p);

    % ---- Optical depth profile (tau integrated from cloud top down) ----
    % tau_lev(i) = optical depth from cloud top to level i
    tau_incr  = (3/2) * (lwc_p ./ re_p) * 1e3 / rho_L .* [0; abs(diff(alt_p))];
    tau_cumTop = cumsum(flipud(tau_incr));   % cumulative from cloud top down
    tau_cumTop = flipud(tau_cumTop);         % restore base-first ordering

    tau_c_total = tau_cumTop(1);             % total cloud optical depth


    % =======================================================================
    % (b) Linear profile retrieval
    % r_e varies linearly from cloud base (r_bot) to cloud top (r_top)
    % =======================================================================
    r_bot_lin = re_p(1);      % cloud base r_e  [µm]
    r_top_lin = re_p(end);    % cloud top  r_e  [µm]

    re_linear_top(kk) = r_top_lin;
    re_linear_bot(kk) = r_bot_lin;

    % Normalised height: 0 at base, 1 at top
    z_norm = (alt_p - alt_p(1)) / max(alt_p(end) - alt_p(1), 1);
    re_profile_lin = r_bot_lin + (r_top_lin - r_bot_lin) .* z_norm(:);

    % LWC consistent with this r_e profile and same total tau
    % tau = (3/2) * integral(LWC / r_e) * (1/rho_L) dz
    % → LWC(z) = r_e(z) * dτ/dz * (2/3) * rho_L
    % We set LWC proportional to r_e (which gives a non-adiabatic shape)
    % scaled so that total tau = tau_c_total
    dz_p     = abs([diff(alt_p(:)); mean(diff(alt_p))]);   % m, layer thickness
    tau_lin_unnorm = (3/2) * (1/rho_L) * 1e3 * sum(dz_p);  % if LWC=r_e everywhere
    % Actually the correct approach: prescribe LWC_lin proportional to r_e_lin
    % such that tau is the same.
    % tau = (3/2)/rho_L * integral(LWC/r_e) dz = (3/2)/rho_L * integral(dz) * N_lev
    % → if LWC ∝ r_e: integral(LWC/r_e) = integral(constant) → tau ∝ cloud depth
    % Better: prescribe a uniform LWC = LWP_target/H (so LWP preserved) and let r_e vary
    % For retrieval comparison we preserve the same tau_c, same LWC distribution (measured)
    lwc_lin = lwc_p(:);   % keep measured LWC; only r_e shape changes

    % Verify that rescaled r_e gives same tau as in-situ
    tau_check_lin = (3/2) * trapz(alt_p, (lwc_lin ./ re_profile_lin) * 1e3 / rho_L);
    % Scale r_e_lin to exactly match tau_c_total
    if tau_check_lin > 0
        scale_r = tau_check_lin / tau_c_total;
        re_profile_lin = re_profile_lin * scale_r;
    end


    % =======================================================================
    % (c) TBLUT-style effective radius
    % r_eff weighted by exp(-tau_above) — mimics 2-band retrieval sensitivity
    % which peaks near cloud top (where most photons are reflected)
    % =======================================================================
    tau_above = tau_cumTop;   % tau from cloud top to each level (0 at top)
    wt_above  = exp(-tau_above);   % exponential weight: largest at cloud top
    re_tblut  = sum(wt_above .* re_p(:) .* tau_incr(:)) / ...
        max(sum(wt_above .* tau_incr(:)), eps);

    re_tblut_all(kk) = re_tblut;

    % Constant r_e = re_tblut, with measured LWC profile
    re_tblut_prof = repmat(re_tblut, N_lev, 1);
    % Scale r_e so that tau matches in-situ tau_c
    tau_check_tblut = (3/2) * trapz(alt_p, (lwc_p ./ re_tblut_prof) * 1e3 / rho_L);
    if tau_check_tblut > 0
        scale_tb = tau_check_tblut / tau_c_total;
        re_tblut_prof = re_tblut_prof * scale_tb;
    end


    % ---- Prepare wc-file arrays (top-to-bottom ordering) ----
    if z_km(2) > z_km(1)
        % currently base-to-top → flip
        z_wc2           = flipud(z_km);
        lwc_wc2         = flipud(lwc_p(:));
        re_wc2_insitu   = flipud(re_p(:));
        re_wc2_linear   = flipud(re_profile_lin);
        re_wc2_tblut    = flipud(re_tblut_prof);
    else
        z_wc2         = z_km;
        lwc_wc2       = lwc_p(:);
        re_wc2_insitu = re_p(:);
        re_wc2_linear = re_profile_lin;
        re_wc2_tblut  = re_tblut_prof;
    end

    if isfield(prof, 'dateOfFlight')
        date_str = prof.dateOfFlight;
    else
        date_str = '2008-11-09';
    end
    if isfield(prof, 'time_utc')
        time_val = mean(prof.time_utc(~isnan(prof.time_utc)), 'omitnan');
    else
        time_val = 120000;
    end

    RT.sza            = sza;
    RT.phi0           = phi0;
    RT.day_of_year    = day_of_year;
    RT.surface_albedo = surface_albedo;
    RT.atm_file       = 'afglus.dat';
    RT.source_file    = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';

    % -- (a) In-situ --
    try
        albedo_retrieval.insitu(kk) = compute_broadband_albedo(...
            re_wc2_insitu, lwc_wc2, z_wc2, ...
            RT, folder_paths, kk, 'retr_insitu', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (in-situ, Part2): %s', kk, ME.message)
    end

    % -- (b) Linear profile retrieval --
    try
        albedo_retrieval.profile(kk) = compute_broadband_albedo(...
            re_wc2_linear, lwc_wc2, z_wc2, ...
            RT, folder_paths, kk, 'retr_linear', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (linear, Part2): %s', kk, ME.message)
    end

    % -- (c) TBLUT --
    try
        albedo_retrieval.tblut(kk) = compute_broadband_albedo(...
            re_wc2_tblut, lwc_wc2, z_wc2, ...
            RT, folder_paths, kk, 'retr_tblut', date_str, time_val, campaign_name);
    catch ME
        warning('Profile %d (TBLUT, Part2): %s', kk, ME.message)
    end

    fprintf('  Part 2: profile %d/%d  |  albedo: insitu=%.4f  linear=%.4f  TBLUT=%.4f\n', ...
        kk, N_use, albedo_retrieval.insitu(kk), ...
        albedo_retrieval.profile(kk), albedo_retrieval.tblut(kk))

end % kk loop Part 2


%% =========================================================================
%  PART 3 — Cloud-top gradient analysis and CTEI
%  =========================================================================
%
%  Compute the CTEI parameter κ for each profile (following
%  oracles_inSitu_analysis.m / Randall 1980, Kuo & Schubert 1988).
%  Then rank profiles by cloud-top r_e gradient and study how:
%    * The gradient magnitude relates to κ (entrainment strength)
%    * The albedo bias (retrieval vs truth) correlates with gradient
%  =========================================================================

fprintf('\n--- Part 3: Cloud-top gradient & CTEI analysis ---\n')

kappa_c = 0.23;          % critical CTEI parameter (Kuo & Schubert 1988)
P0      = 1e5;           % reference pressure [Pa]

kappa_CTEI  = NaN(1, N_use);
DeltaT_inv  = NaN(1, N_use);

for kk = 1:N_use

    nn   = idx_use(kk);
    prof = ensemble_profiles{nn};

    [~, lwc_p, ~, alt_p, T_K_p, P_Pa_p, ~] = orientProfile(prof);

    valid_lev = (lwc_p > 0) & ~isnan(T_K_p) & ~isnan(P_Pa_p);
    if sum(valid_lev) < 4
        continue
    end
    lwc_p  = lwc_p(valid_lev);
    alt_p  = alt_p(valid_lev);
    T_K_p  = T_K_p(valid_lev);
    P_Pa_p = P_Pa_p(valid_lev);

    T_base  = T_K_p(1);
    P_base  = P_Pa_p(1);
    T_top   = T_K_p(end);
    P_top   = P_Pa_p(end);

    % Saturation vapour pressure [Pa] via Tetens
    e_s_base = 611.2 * exp(17.67 * (T_base - 273.15) / (T_base - 29.65));
    r_s_base = epsilon * e_s_base / max(P_base - e_s_base, 1);

    % Moist adiabatic lapse rate at cloud base [K/m]
    SALR = g * (1 + L_v * r_s_base / (R_d * T_base)) / ...
        (c_pd + L_v^2 * r_s_base / (R_v * T_base^2));

    H_total = abs(alt_p(end) - alt_p(1));
    T_moist_adiabat_top = T_base - SALR * H_total;

    DeltaT_inv(kk) = T_top - T_moist_adiabat_top;   % K; > 0 → warmer than expected

    if DeltaT_inv(kk) > 0
        e_s_top = 611.2 * exp(17.67 * (T_top - 273.15) / (T_top - 29.65));
        q_s_top = epsilon * e_s_top / max(P_top - e_s_top, 1);

        T_FT_est  = T_top + DeltaT_inv(kk);
        e_s_FT    = 611.2 * exp(17.67 * (T_FT_est - 273.15) / (T_FT_est - 29.65));
        q_s_FT    = epsilon * e_s_FT / max(P_top - e_s_FT, 1);

        Delta_q_s = q_s_top - q_s_FT;   % < 0: free troposphere is drier

        if abs(Delta_q_s) > 1e-6
            kappa_CTEI(kk) = c_pd * DeltaT_inv(kk) / (L_v * abs(Delta_q_s));
        end
    end

end






%% =========================================================================
%  PART 4 — Fixed-tau_c albedo comparison - different from Borg and Bennartz
%  =========================================================================
%
%  For each profile the same measured cloud optical depth is
%  kept.  Three effective-radius profile shapes are imposed:
%    (a) in-situ measured r_e
%    (b) adiabatic shape  r_e ∝ (z – z_base)^(1/3), scaled to same LWP-tau
%    (c) constant r_e    = LWC-weighted mean of in-situ r_e
%  =========================================================================



%% =========================================================================
%  SAVE RESULTS
%  =========================================================================

results.campaign_name     = campaign_name;
results.sza               = sza;
results.day_of_year       = day_of_year;
results.idx_use           = idx_use;
results.N_use             = N_use;
results.LWP_all           = LWP_all;
results.cloud_depth_all   = cloud_depth_all;
results.cloudTop_re_grad  = cloudTop_re_grad;
results.re_top_insitu     = re_top_insitu;
results.re_base_insitu    = re_base_insitu;
results.re_tblut_all      = re_tblut_all;
results.re_linear_top     = re_linear_top;
results.re_linear_bot     = re_linear_bot;
results.tau_insitu        = tau_insitu;
results.tau_adiabatic     = tau_adiabatic;
results.tau_constant      = tau_constant;
results.albedo_fixedLWP   = albedo_fixedLWP;
results.albedo_retrieval  = albedo_retrieval;
results.kappa_CTEI        = kappa_CTEI;
results.DeltaT_inv        = DeltaT_inv;

save_filename = [foldername_save, 'cloud_albedo_profile_analysis_', ...
    campaign_name, '_sza', num2str(sza), '_', datestr(now, 'dd-mmm-yyyy'), '.mat'];
save(save_filename, 'results')
fprintf('\nResults saved to:\n  %s\n', save_filename)


%% =========================================================================
%  FIGURES
%  =========================================================================

% ---- Color palette ----
clr_insitu    = [0.12, 0.47, 0.71];   % blue
clr_adiabatic = [0.20, 0.63, 0.17];   % green
clr_constant  = [0.89, 0.10, 0.11];   % red
clr_tblut     = [1.00, 0.50, 0.05];   % orange
clr_linear    = [0.54, 0.17, 0.89];   % purple

% ---- Mask out NaN rows ----
ok1 = ~isnan(albedo_fixedLWP.insitu) & ~isnan(albedo_fixedLWP.adiabatic) & ...
    ~isnan(albedo_fixedLWP.constant);
ok2 = ~isnan(albedo_retrieval.insitu) & ~isnan(albedo_retrieval.profile) & ...
    ~isnan(albedo_retrieval.tblut);
ok3 = ok1 & ~isnan(kappa_CTEI) & ~isnan(cloudTop_re_grad);


% -----------------------------------------------------------------------
% Figure 1 — Fixed LWP: scatter of albedos for each profile type vs LWP
% -----------------------------------------------------------------------
figure('Position', [50 50 1300 450])
subplot(1,3,1)
scatter(LWP_all(ok1), albedo_fixedLWP.insitu(ok1), 30, clr_insitu, 'filled', ...
    'MarkerFaceAlpha', 0.7)
xlabel('LWP (g m^{-2})'); ylabel('Broadband albedo')
title('In-situ r_e profile')
grid on; grid minor; ylim([0 1])

subplot(1,3,2)
scatter(LWP_all(ok1), albedo_fixedLWP.adiabatic(ok1), 30, clr_adiabatic, 'filled', ...
    'MarkerFaceAlpha', 0.7)
xlabel('LWP (g m^{-2})');
title('Adiabatic r_e profile')
grid on; grid minor; ylim([0 1])

subplot(1,3,3)
scatter(LWP_all(ok1), albedo_fixedLWP.constant(ok1), 30, clr_constant, 'filled', ...
    'MarkerFaceAlpha', 0.7)
xlabel('LWP (g m^{-2})');
title('Constant r_e profile')
grid on; grid minor; ylim([0 1])

sgtitle(['Fixed LWP broadband albedo — ', campaign_name, ...
    '   SZA = ', num2str(sza), '\circ'], 'FontSize', 13)
set(gcf, 'Name', 'fig1_fixedLWP_albedo_vs_LWP')


% -----------------------------------------------------------------------
% Figure 2 — Fixed LWP: albedo differences between profile types
% -----------------------------------------------------------------------
da_ad_vs_in  = albedo_fixedLWP.adiabatic - albedo_fixedLWP.insitu;
da_const_vs_in = albedo_fixedLWP.constant - albedo_fixedLWP.insitu;

figure('Position', [50 550 900 500])
hold on
histogram(da_ad_vs_in(ok1), 'BinWidth', 0.005, ...
    'FaceColor', clr_adiabatic, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'Adiabatic – In-situ')
histogram(da_const_vs_in(ok1), 'BinWidth', 0.005, ...
    'FaceColor', clr_constant, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'Constant – In-situ')
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xlabel('\Delta Broadband Albedo (relative to in-situ)')
ylabel('Relative frequency')
title(['Fixed LWP: albedo difference — ', campaign_name])
legend('Location', 'best')
grid on; grid minor
set(gcf, 'Name', 'fig2_fixedLWP_albedo_diff_hist')


% -----------------------------------------------------------------------
% Figure 3 — Fixed LWP: albedo vs albedo scatter (truth vs idealized)
% -----------------------------------------------------------------------
figure('Position', [50 50 900 420])
subplot(1,2,1)
hold on
plot([0 1], [0 1], 'k--', 'LineWidth', 1)
scatter(albedo_fixedLWP.insitu(ok1), albedo_fixedLWP.adiabatic(ok1), ...
    30, clr_adiabatic, 'filled', 'MarkerFaceAlpha', 0.7)
xlabel('Albedo — in-situ'); ylabel('Albedo — adiabatic')
title('Adiabatic vs. In-situ (fixed LWP)')
axis equal; xlim([0 1]); ylim([0 1])
grid on; grid minor

subplot(1,2,2)
hold on
plot([0 1], [0 1], 'k--', 'LineWidth', 1)
scatter(albedo_fixedLWP.insitu(ok1), albedo_fixedLWP.constant(ok1), ...
    30, clr_constant, 'filled', 'MarkerFaceAlpha', 0.7)
xlabel('Albedo — in-situ'); ylabel('Albedo — constant r_e')
title('Constant vs. In-situ (fixed LWP)')
axis equal; xlim([0 1]); ylim([0 1])
grid on; grid minor

sgtitle(['Fixed LWP: albedo scatter — ', campaign_name], 'FontSize', 12)
set(gcf, 'Name', 'fig3_fixedLWP_albedo_scatter')


% -----------------------------------------------------------------------
% Figure 4 — Retrieval comparison: albedo differences vs truth
% -----------------------------------------------------------------------
da_lin_vs_in  = albedo_retrieval.profile - albedo_retrieval.insitu;
da_tblut_vs_in = albedo_retrieval.tblut  - albedo_retrieval.insitu;

figure('Position', [50 550 900 500])
hold on
histogram(da_lin_vs_in(ok2), 'BinWidth', 0.005, ...
    'FaceColor', clr_linear, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'Linear profile – In-situ')
histogram(da_tblut_vs_in(ok2), 'BinWidth', 0.005, ...
    'FaceColor', clr_tblut, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'TBLUT – In-situ')
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xlabel('\Delta Broadband Albedo (relative to true in-situ)')
ylabel('Relative frequency')
title(['Retrieval albedo bias — ', campaign_name, '   SZA = ', num2str(sza), '\circ'])
legend('Location', 'best')
grid on; grid minor
set(gcf, 'Name', 'fig4_retrieval_albedo_diff_hist')


% -----------------------------------------------------------------------
% Figure 5 — Cloud-top gradient: scatter vs kappa_CTEI
% -----------------------------------------------------------------------
figure('Position', [50 50 700 550])
hold on
scatter(cloudTop_re_grad(ok3), kappa_CTEI(ok3), 50, ...
    LWP_all(ok3), 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
cb = colorbar;
cb.Label.String = 'LWP (g m^{-2})';
yline(kappa_c, 'r--', 'LineWidth', 1.5, 'DisplayName', ['\kappa_c = ', num2str(kappa_c)])
xlabel('Cloud-top r_e gradient  dr_e/dz  (μm m^{-1})')
ylabel('CTEI parameter  \kappa')
title(['Cloud-top gradient vs. CTEI — ', campaign_name])
legend('Location', 'best')
grid on; grid minor
set(gcf, 'Name', 'fig5_cloudTop_gradient_vs_CTEI')


% -----------------------------------------------------------------------
% Figure 6 — Albedo bias (TBLUT – in-situ) vs cloud-top r_e gradient
% -----------------------------------------------------------------------
ok6 = ok2 & ~isnan(cloudTop_re_grad);

figure('Position', [50 600 700 500])
hold on
scatter(cloudTop_re_grad(ok6), da_tblut_vs_in(ok6), 50, ...
    kappa_CTEI(ok6), 'filled', 'MarkerFaceAlpha', 0.8)
cb = colorbar;
cb.Label.String = '\kappa_{CTEI}';
yline(0, 'k--', 'LineWidth', 1.2)
xlabel('Cloud-top r_e gradient  dr_e/dz  (μm m^{-1})')
ylabel('\Delta Albedo = TBLUT – In-situ')
title(['TBLUT albedo bias vs. cloud-top gradient — ', campaign_name])
grid on; grid minor
set(gcf, 'Name', 'fig6_TBLUT_bias_vs_cloudTopGrad')


% -----------------------------------------------------------------------
% Figure 7 — Summary: all 3 retrieval albedos vs in-situ for best/worst cases
%
%  Identify the profiles with the largest cloud-top r_e gradient
%  (positive: r_e decreasing toward top, negative: r_e increasing toward top)
%  and show the profile shapes + albedo comparison for those profiles.
% -----------------------------------------------------------------------
[~, idx_sorted] = sort(abs(cloudTop_re_grad), 'descend', 'MissingPlacement', 'last');
N_highlight = min(6, sum(~isnan(cloudTop_re_grad)));

figure('Position', [50 50 1300 650])

for jj = 1:N_highlight
    kk = idx_sorted(jj);
    if ~ok2(kk) || isnan(cloudTop_re_grad(kk))
        continue
    end
    nn   = idx_use(kk);
    prof = ensemble_profiles{nn};

    [re_p, lwc_p, ~, alt_p, ~, ~, ~] = orientProfile(prof);
    valid_lev = (re_p > 0.5) & (lwc_p > 0) & ~isnan(re_p);
    re_p  = re_p(valid_lev);
    alt_p = alt_p(valid_lev);

    z_norm_prof = (alt_p - alt_p(1)) / max(alt_p(end) - alt_p(1), 1);

    subplot(2, N_highlight, jj)
    plot(re_p, z_norm_prof, 'Color', clr_insitu, 'LineWidth', 1.5)
    xlabel('r_e (µm)'); ylabel('Normalised height')
    title(sprintf('Prof %d  |  grad = %.3f µm/m', kk, cloudTop_re_grad(kk)), ...
        'FontSize', 8)
    grid on

    subplot(2, N_highlight, N_highlight + jj)
    bar([albedo_retrieval.insitu(kk), ...
        albedo_retrieval.profile(kk), ...
        albedo_retrieval.tblut(kk)], 0.5, 'FaceColor', 'flat', ...
        'CData', [clr_insitu; clr_linear; clr_tblut])
    set(gca, 'XTickLabel', {'In-situ', 'Linear', 'TBLUT'}, 'XTickLabelRotation', 30)
    ylabel('Broadband albedo')
    ylim([0 1]); grid on
end

sgtitle(['Profiles with largest cloud-top r_e gradients — ', campaign_name], ...
    'FontSize', 12)
set(gcf, 'Name', 'fig7_highlight_profiles')


fprintf('\n=== Analysis complete. Results saved. ===\n')


%% =========================================================================
%  LOCAL HELPER: orientProfile
%  (Duplicated here so the script is self-contained if the function is not
%  on the path.  If you already have orientProfile.m on your path, this
%  local copy is shadowed by it.)
%  =========================================================================

function [re, lwc, Nc, alt, T_K, P_Pa, rv] = orientProfile(prof)
% ORIENTPROFILE  Orient a vertical profile to cloud-base-first order.
dz_dt = mean(diff(prof.altitude) ./ diff(prof.time));

if dz_dt > 0

    if isfield(prof, 're')
        re   = prof.re;
    elseif isfield(prof, 're_CDP')
        re = prof.re_CDP;
    end

    lwc  = prof.lwc;
    Nc   = prof.total_Nc;
    alt  = prof.altitude;
    T_K  = prof.temp + 273.15;
    P_Pa = prof.pres * 100;
    % rv   = prof.rv;
    rv = [];
else

    if isfield(prof, 're')
        re   = fliplr(prof.re);
    elseif isfield(prof, 're_CDP')
        re = fliplr(prof.re_CDP);
    end

    lwc  = fliplr(prof.lwc);
    Nc   = fliplr(prof.total_Nc);
    alt  = fliplr(prof.altitude);
    T_K  = fliplr(prof.temp + 273.15);
    P_Pa = fliplr(prof.pres * 100);
    % rv   = fliplr(prof.rv);
    rv = [];
end

re   = re(:)';
lwc  = lwc(:)';
Nc   = Nc(:)';
alt  = alt(:)';
T_K  = T_K(:)';
P_Pa = P_Pa(:)';
% rv   = rv(:)';
rv = [];
end
