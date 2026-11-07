%% Estimate the adiabaticity of one in-situ vertical cloud profile (Wood 2005)
%
% DESCRIPTION:
%   Two adiabaticity estimates for a single aircraft vertical profile through
%   a liquid cloud, following Wood (2005, J. Atmos. Sci., 62, 3011-3033,
%   section 3b, p. 3019):
%
%   (1) PROFILE method. The adiabatic liquid water content at height z is
%
%           q_ad(z) = Gamma_ad * (z - z_cb)
%
%       and the per-level adiabaticity is f_ad(z) = LWC(z) / q_ad(z)
%       (Wood's Fig. 4c). A scalar summary is the mean of f_ad(z) over the
%       levels with normalized height z* = (z - z_cb)/h >= zstar_min; the
%       levels nearest cloud base are excluded because q_ad -> 0 there and
%       the ratio is dominated by noise.
%
%   (2) PATH method. The cloud's overall adiabaticity is
%
%           f_ad = LWP_obs / LWP_ad,     LWP_ad = (1/2) * Gamma_ad * h^2
%
%       with h = z_top - z_cb (Wood's Table 3).
%
%   If LWC(z) = f * q_ad(z) with a constant f, both methods return exactly f.
%   They differ when the per-level ratio varies with height: the path method
%   weights the upper levels more heavily because q_ad is largest there, so
%   an entrainment-depleted cloud top pulls f_ad_lwp below f_ad_profile.
%
%   Gamma_ad(T, p) is the adiabatic rate of increase of LWC with height.
%   Wood: "A single value of Gamma_ad is used in each case, calculated using
%   the mean in-cloud temperature and pressure." That is done here: T and p
%   are averaged over the sampled cloud depth (trapezoid weighting in altitude)
%   and Gamma_ad is evaluated once at that mean state.
%
% CLOUD BASE (base_method = 'extrapolate', the default):
%   The profile finders keep only samples that pass their LWC and Nc thresholds,
%   so the lowest kept sample sits some distance above the true cloud base, and
%   right at cloud base the measurements are ambiguous because the droplet number
%   density is small and the probes are noisy. Both methods above depend on z_cb,
%   and f_ad_lwp scales as 1/h^2, so the base is estimated from the cloud itself:
%
%   The adiabatic model is linear in height above cloud base,
%       w_l = Gamma_l (z - z_B)          Albrecht et al. (1990, GRL 17, 89-92,
%                                        doi:10.1029/GL017i001p00089)
%   and a sub-adiabatic cloud with a height-independent f_ad is linear too, with
%   slope f_ad * Gamma_ad. So a straight line is fit to LWC(z) over the bottom
%   fit_fraction of the sampled cloud depth, and its z-intercept (LWC = 0) is
%   taken as z_cb. The line is extrapolated below the lowest kept sample; it
%   is not fit to anything there.
%
%   FALLBACK. If the fitted slope is not physical, the base is placed one
%   adiabatic distance below the lowest sample, assuming that sample sits at
%   the finder's LWC threshold:
%       z_cb = z_lowest - lwc_threshold / Gamma_ad
%   The slope is treated as not physical when
%       slope <= slope_limits_over_gamma_ad(1) * Gamma_ad   (no usable upward trend)
%       slope >  slope_limits_over_gamma_ad(2) * Gamma_ad   (far steeper than adiabatic)
%   or when fewer than 3 samples are available to fit. (With LWC >= 0 and a
%   positive slope, the least-squares line passes through the centroid of the
%   fitted samples, so it always reaches zero below the top of the fitted layer;
%   no separate check is needed for that.)
%
%   WHY A LOWER LIMIT ABOVE ZERO. The intercept distance is LWC_fit(z_lowest) /
%   slope, so a slope that is positive but tiny sends the base very far down.
%   With only "slope <= 0" excluded, the bottom-third fits put the base as much
%   as 1.2 km (VOCALS-REx) and 7 km (ORACLES) below the cloud. Requiring at least
%   0.25 Gamma_ad (a quarter of adiabatic) limits the extrapolation to 161 m
%   (VOCALS-REx) and 292 m (ORACLES) and changes the median f_ad by 0.02 or
%   less in either campaign. The upper limit of 3 Gamma_ad rejects a line that
%   rises faster than any adiabatic parcel could.
%
%   For reference, the adiabatic rate of increase of LWC with height is usually
%   2.0-2.8 g/m^3 per km in marine stratus and stratocumulus, and 2.1-2.6 (median
%   2.24) over the southeast Pacific during EPIC 2001 (Zhou et al. 2006, GRL 33,
%   L02807, doi:10.1029/2005GL024796, para. 8, citing Albrecht et al. 1990,
%   Austin et al. 1995 and Wood and Taylor 2001). The limits above scale with each
%   profile's own Gamma_ad, which depends on its mean (T, p); the colder VOCALS-REx
%   clouds (mean in-cloud T 3.6-11.6 deg C) give 1.72-2.21, partly below that range.
%
%   THE BASE IS NEVER ABOVE THE LOWEST SAMPLE:
%       z_cb = min(line intercept, z_lowest)
%   A profile can end in a few samples with almost no cloud water (a descent,
%   for example, continues until LWC and Nc are both below threshold), and the
%   line then reaches zero inside the profile. In that case the lowest sample is
%   the base and every sample counts toward LWP_obs (A. Buggee, Nov 2026). On
%   the 14 such profiles in VOCALS-REx and ORACLES this lowers f_ad by a median
%   factor of 0.98 (0.93-1.00) relative to using the intercept.
%   base_source then reads 'lowest sample, intercept above it'.
%
%   LWP_obs IS INTEGRATED OVER THE SAME LAYER AS LWP_ad. When the base is below
%   the lowest sample, the unsampled wedge between them is added using the same
%   line: (1/2) LWC_line(z_lowest) (z_lowest - z_cb), or (1/2) lwc_threshold dz
%   for the fallback. The wedge is an estimate, not a measurement; it is a
%   median of 1-2% of LWP_obs.
%
%   base_method = 'lowest_sample' takes the lowest sample as the base with no
%   correction. The VOCALS-REx finder also stores a measured base for comparison
%   (find_insitu_cloud_base.m, called with 'lowest_sample' on an extended profile).
%
% WHICH LWC. Everything above uses lwc_g_m3, normally the total LWC (cloud
%   droplets plus drizzle drops), the same quantity as the profile's .lwp field.
%   If the droplet-mode LWC is also supplied, the *_droplet fields repeat the
%   numerator with drizzle drops left out; the cloud base is not changed by it.
%
% GAMMA_AD. Wood does not give a formula for Gamma_ad. It follows from
%   conservation of moist static energy along a saturated adiabat,
%
%       c_p dT + g dz + L_v dq* = 0   =>   dq_l/dz = (c_p / L_v) (Gamma_d - Gamma_m)
%
%   so that, per unit volume of air,
%
%       Gamma_ad = rho_air * (c_p / L_v) * (Gamma_d - Gamma_m)
%
%   where Gamma_m is the moist adiabatic lapse rate of Randall (2009, "The
%   Moist Adiabatic Lapse Rate", QuickStudies for Graduate Students in
%   Atmospheric Science, eq. 36). Gamma_m is a state function of the local
%   (T, p) -- the slope of the moist adiabat through that point -- NOT a
%   finite difference of the observed temperature across the cloud. The
%   same expression appears in Albrecht et al. (1990) and Brenguier (1991,
%   JAS 48, 264).
%
%   The formulas, the constants and the cloud-base rule are identical to the
%   Python implementation in
%       lasp-CU-paper-3/plot_per_cloud_rmse_kfold.py
%   (saturation_vapor_pressure_Pa, saturation_specific_humidity,
%   moist_adiabatic_lapse_rate, gamma_ad_kg_per_m4, extrapolated_cloud_base,
%   wood_adiabaticity) so the in-situ and synthetic pipelines are directly
%   comparable. If one is changed, change the other.
%
% INPUTS:
%   altitude_m - altitude of each sample [m], in sampling order (a profile
%                may be flown climbing or descending; either is fine)
%   lwc_g_m3   - liquid water content of each sample [g/m^3]
%   temp_C     - ambient air temperature of each sample [deg C]
%   pres_hPa   - ambient static pressure of each sample [hPa]
%
%   Optional name-value arguments:
%   'base_method'                 'extrapolate' (default) or 'lowest_sample'
%   'lwc_threshold_g_m3'          the profile finder's LWC threshold, used by the
%                                 fallback [g/m^3] (default 0.05)
%   'fit_fraction'                bottom fraction of the sampled depth fit by the
%                                 line [-] (default 1/3)
%   'slope_limits_over_gamma_ad'  [lower, upper] limits on the fitted slope as
%                                 multiples of Gamma_ad [-] (default [0.25, 3]).
%                                 [0, Inf] keeps only the "slope <= 0" test.
%   'zstar_min'                   lower bound on z* when averaging f_ad(z) into
%                                 f_ad_profile [-] (default 0.1)
%   'lwc_droplet_g_m3'            droplet-mode LWC of each sample [g/m^3]; fills
%                                 the *_droplet fields (default: not supplied)
%
% OUTPUT:
%   adiab - structure with fields
%       .T_mean_C              mean in-cloud temperature used for Gamma_ad [deg C]
%       .P_mean_hPa            mean in-cloud pressure used for Gamma_ad [hPa]
%       .Gamma_m_K_per_km      moist adiabatic lapse rate at (T_mean, P_mean) [K/km]
%       .Gamma_ad_kg_per_m4    adiabatic LWC lapse rate [kg/m^3 per m]
%       .Gamma_ad_g_m3_per_km  same, in the units people quote [g/m^3 per km]
%       .base_method           'extrapolate' or 'lowest_sample'
%       .base_source           'fit intercept', 'fallback', 'lowest sample, intercept
%                              above it', or 'lowest sample' (base_method 'lowest_sample')
%       .n_fit_samples         samples in the bottom fit_fraction of the cloud
%       .fit_slope_g_m3_per_km slope of the line fit to LWC(z) [g/m^3 per km]
%       .fit_slope_over_gamma_ad   that slope as a multiple of Gamma_ad [-]
%       .lwc_fit_at_lowest_g_m3    the line evaluated at the lowest sample [g/m^3]
%       .z_lowest_sample_m     altitude of the lowest sample [m]
%       .dz_base_m             z_lowest_sample - z_cb [m]; >= 0, and 0 when the
%                              lowest sample is the base
%       .z_cb_m                cloud base used [m]
%       .z_top_m               cloud top = highest sample altitude [m]
%       .h_m                   cloud thickness z_top - z_cb [m]
%       .z_star                (z - z_cb)/h for each sample, sampling order [-]
%       .q_ad_g_m3             adiabatic LWC at each sample, sampling order [g/m^3]
%       .f_ad_z                LWC(z)/q_ad(z) at each sample; NaN where q_ad <= 0 [-]
%       .f_ad_profile          method (1) scalar: mean of f_ad_z over z* >= zstar_min [-]
%       .f_ad_z_droplet        droplet-mode LWC(z)/q_ad(z) at each sample [-]
%       .f_ad_profile_droplet  mean of f_ad_z_droplet over z* >= zstar_min [-]
%       .lwp_measured_g_m2     |trapz(altitude, LWC)| over every sample, the same
%                              integral as the profile's .lwp field [g/m^2]
%       .lwp_wedge_g_m2        unsampled layer between z_cb and the lowest sample [g/m^2]
%       .lwp_obs_g_m2          LWP from z_cb to cloud top, wedge included [g/m^2]
%       .lwp_ad_g_m2           (1/2) Gamma_ad h^2 [g/m^2]
%       .f_ad_lwp              method (2) scalar: lwp_obs / lwp_ad [-]
%       .lwp_droplet_g_m2      droplet-mode LWP over the same layer [g/m^2]
%       .f_ad_lwp_droplet      lwp_droplet / lwp_ad [-]
%       .h_lowest_sample_m     z_top - lowest sample [m]
%       .lwp_ad_lowest_sample_g_m2   (1/2) Gamma_ad h_lowest_sample^2 [g/m^2]
%       .f_ad_lwp_lowest_sample      lwp_measured / lwp_ad_lowest_sample [-]
%       .lwc_threshold_g_m3, .fit_fraction, .slope_limits_over_gamma_ad,
%       .zstar_min             settings used
%
%   Every numeric field is NaN when the profile has no vertical extent, has
%   fewer than two valid (T, p) samples, or yields Gamma_ad <= 0.
%
% SEE ALSO:
%   find_verticalProfiles_VOCALS_REx_ver2, find_verticalProfiles_ORACLES,
%   find_insitu_cloud_base, plot_ensemble_adiabaticity,
%   plot_adiabaticity_base_comparison
%
% By Andrew John Buggee
%%

function adiab = compute_adiabaticity_Wood2005(altitude_m, lwc_g_m3, temp_C, pres_hPa, options)

arguments
    altitude_m double
    lwc_g_m3 double
    temp_C double
    pres_hPa double
    options.base_method (1,:) char {mustBeMember(options.base_method, {'extrapolate', 'lowest_sample'})} = 'extrapolate'
    options.lwc_threshold_g_m3 (1,1) double {mustBePositive} = 0.05
    options.fit_fraction (1,1) double {mustBeInRange(options.fit_fraction, 0, 1, 'exclude-lower')} = 1/3
    options.slope_limits_over_gamma_ad (1,2) double {mustBeNonnegative} = [0.25, 3]
    options.zstar_min (1,1) double = 0.1
    options.lwc_droplet_g_m3 double = []
end

base_method = options.base_method;
zstar_min   = options.zstar_min;

% Work in double row vectors regardless of how the profile was stored
z_m   = altitude_m(:)';                  % m
lwc   = lwc_g_m3(:)';                    % g/m^3
T_K   = temp_C(:)' + 273.15;             % K
P_Pa  = pres_hPa(:)' .* 100;             % Pa

if numel(z_m) ~= numel(lwc) || numel(z_m) ~= numel(T_K) || numel(z_m) ~= numel(P_Pa)
    error([newline, 'altitude, lwc, temp and pres must have the same number of samples.', newline])
end

% The droplet-mode LWC is optional. Without it the *_droplet fields stay NaN.
has_droplet_lwc = ~isempty(options.lwc_droplet_g_m3);

if has_droplet_lwc
    lwc_droplet = options.lwc_droplet_g_m3(:)';      % g/m^3

    if numel(lwc_droplet) ~= numel(z_m)
        error([newline, 'lwc_droplet_g_m3 must have the same number of samples as altitude.', newline])
    end
end

adiab = nan_output(z_m, options);

% ------------------------------------------------------------------------
% Sampled geometry. The lowest and highest samples span the .cloud_depth
% field of the profile structures. The cloud base used for the
% adiabaticity is set further down, once Gamma_ad is known.
% ------------------------------------------------------------------------
z_lowest_m = min(z_m);
z_top_m    = max(z_m);
h_sampled_m = z_top_m - z_lowest_m;

adiab.z_lowest_sample_m = z_lowest_m;
adiab.z_top_m           = z_top_m;
adiab.h_lowest_sample_m = h_sampled_m;

if ~(h_sampled_m > 0)
    return
end

% ------------------------------------------------------------------------
% Mean in-cloud temperature and pressure (Wood 2005). Trapezoid weighting
% in altitude over the sampled profile, divided by the net altitude change,
% which is the depth-weighted mean even if the aircraft altitude wobbles.
% Samples where T or p are missing (netCDF fill -> NaN) are dropped.
% ------------------------------------------------------------------------
valid = isfinite(T_K) & isfinite(P_Pa) & isfinite(z_m);

if sum(valid) < 2 || abs(z_m(find(valid, 1, 'last')) - z_m(find(valid, 1, 'first'))) <= 0
    return
end

dz_net = z_m(find(valid, 1, 'last')) - z_m(find(valid, 1, 'first'));       % m, signed
T_mean_K  = trapz(z_m(valid), T_K(valid))  ./ dz_net;                       % K
P_mean_Pa = trapz(z_m(valid), P_Pa(valid)) ./ dz_net;                       % Pa

adiab.T_mean_C   = T_mean_K - 273.15;
adiab.P_mean_hPa = P_mean_Pa ./ 100;

% ------------------------------------------------------------------------
% Adiabatic LWC lapse rate at the mean in-cloud state
% ------------------------------------------------------------------------
Gamma_m   = moist_adiabatic_lapse_rate(T_mean_K, P_mean_Pa);              % K/m
Gamma_ad  = gamma_ad_kg_per_m4(T_mean_K, P_mean_Pa);                       % kg/m^3 per m

adiab.Gamma_m_K_per_km     = Gamma_m .* 1e3;
adiab.Gamma_ad_kg_per_m4   = Gamma_ad;
adiab.Gamma_ad_g_m3_per_km = Gamma_ad .* 1e6;      % kg/m^4 -> g/m^3 per km

if ~(Gamma_ad > 0)
    return
end

Gamma_ad_g_m3_per_m = Gamma_ad .* 1e3;             % kg/m^4 -> g/m^3 per m

% ------------------------------------------------------------------------
% Straight-line fit to LWC(z) over the bottom of the sampled cloud. Heights are
% measured from the lowest sample, so the fitted intercept term is the line's
% LWC at the lowest sample. Computed for both base methods, as a diagnostic.
% The 1e-4 tolerance (3 cm in a 300-m cloud) keeps a level that sits exactly on
% the boundary: evenly spaced levels stored as float32, like the 7-level synthetic
% profiles, otherwise land ~1e-6 above it and drop out of the fit.
% ------------------------------------------------------------------------
lwp_valid = isfinite(lwc) & isfinite(z_m);

z_star_sampled = (z_m - z_lowest_m) ./ h_sampled_m;
in_fit = lwp_valid & z_star_sampled <= options.fit_fraction + 1e-4;

adiab.n_fit_samples = sum(in_fit);

if sum(in_fit) >= 3 && std(z_m(in_fit)) > 0
    p_fit = polyfit(z_m(in_fit) - z_lowest_m, lwc(in_fit), 1);
    slope_g_m3_per_m       = p_fit(1);              % g/m^3 per m
    lwc_line_at_lowest_g_m3 = p_fit(2);             % g/m^3
else
    slope_g_m3_per_m       = NaN;
    lwc_line_at_lowest_g_m3 = NaN;
end

slope_over_gamma_ad = slope_g_m3_per_m ./ Gamma_ad_g_m3_per_m;

adiab.fit_slope_g_m3_per_km   = slope_g_m3_per_m .* 1e3;
adiab.fit_slope_over_gamma_ad = slope_over_gamma_ad;
adiab.lwc_fit_at_lowest_g_m3  = lwc_line_at_lowest_g_m3;

% ------------------------------------------------------------------------
% Cloud base
% ------------------------------------------------------------------------
switch base_method

    case 'lowest_sample'

        dz_base_m = 0;
        lwc_at_lowest_for_wedge = 0;
        adiab.base_source = 'lowest sample';

    case 'extrapolate'

        % Is the fitted slope physical? (NaN fails every comparison, so a profile
        % with too few samples to fit also lands on the fallback.)
        slope_is_physical = slope_over_gamma_ad > options.slope_limits_over_gamma_ad(1) && ...
                            slope_over_gamma_ad <= options.slope_limits_over_gamma_ad(2);

        if slope_is_physical

            % Distance from the lowest sample down to where the line reaches LWC = 0.
            % Negative when the line reaches zero inside the profile.
            dz_base_m = lwc_line_at_lowest_g_m3 ./ slope_g_m3_per_m;                    % m

            if dz_base_m >= 0
                lwc_at_lowest_for_wedge = lwc_line_at_lowest_g_m3;                       % g/m^3
                adiab.base_source = 'fit intercept';
            else
                % The line reaches LWC = 0 above the lowest sample. The base is never
                % put above a measured sample, so the lowest sample is the base.
                dz_base_m = 0;
                lwc_at_lowest_for_wedge = 0;
                adiab.base_source = 'lowest sample, intercept above it';
            end

        else
            % Fallback: the lowest sample is taken to sit at the finder's LWC
            % threshold, with an adiabatic layer below it
            dz_base_m = options.lwc_threshold_g_m3 ./ Gamma_ad_g_m3_per_m;              % m
            lwc_at_lowest_for_wedge = options.lwc_threshold_g_m3;                       % g/m^3
            adiab.base_source = 'fallback';
        end

end

z_cb_m = z_lowest_m - dz_base_m;
h_m    = z_top_m - z_cb_m;

adiab.dz_base_m = dz_base_m;
adiab.z_cb_m    = z_cb_m;
adiab.h_m       = h_m;
adiab.z_star    = (z_m - z_cb_m) ./ h_m;

% ------------------------------------------------------------------------
% Method (1): q_ad(z) = Gamma_ad (z - z_cb) and the per-level ratio
% ------------------------------------------------------------------------
q_ad = Gamma_ad_g_m3_per_m .* (z_m - z_cb_m);      % g/m^3, sampling order
f_ad_z = lwc ./ q_ad;
f_ad_z(~(q_ad > 0)) = NaN;                         % at or below cloud base

adiab.q_ad_g_m3 = q_ad;
adiab.f_ad_z    = f_ad_z;

use_for_mean = adiab.z_star >= zstar_min & isfinite(f_ad_z);
if any(use_for_mean)
    adiab.f_ad_profile = mean(f_ad_z(use_for_mean));
end

% Same ratio using the droplet-mode LWC only
if has_droplet_lwc

    f_ad_z_droplet = lwc_droplet ./ q_ad;
    f_ad_z_droplet(~(q_ad > 0)) = NaN;
    adiab.f_ad_z_droplet = f_ad_z_droplet;

    use_for_mean = adiab.z_star >= zstar_min & isfinite(f_ad_z_droplet);
    if any(use_for_mean)
        adiab.f_ad_profile_droplet = mean(f_ad_z_droplet(use_for_mean));
    end

end

% ------------------------------------------------------------------------
% Method (2): f_ad = LWP_obs / LWP_ad, both over the layer from z_cb to cloud top.
% The base is never above the lowest sample, so every sample is in that layer and
% lwp_measured is the same integral as the profile's .lwp field. When z_cb is
% below the lowest sample, the unsampled wedge is added: the line (or the
% fallback's adiabat) is triangular there, so its integral is
% (1/2) LWC(z_lowest) (z_lowest - z_cb).
% ------------------------------------------------------------------------
lwp_measured = abs(trapz(z_m(lwp_valid), lwc(lwp_valid)));                 % g/m^2
lwp_wedge    = 0.5 .* lwc_at_lowest_for_wedge .* dz_base_m;                % g/m^2

adiab.lwp_measured_g_m2 = lwp_measured;
adiab.lwp_wedge_g_m2    = lwp_wedge;
adiab.lwp_obs_g_m2      = lwp_measured + lwp_wedge;
adiab.lwp_ad_g_m2       = 0.5 .* Gamma_ad_g_m3_per_m .* h_m.^2;             % g/m^2
adiab.f_ad_lwp          = adiab.lwp_obs_g_m2 ./ adiab.lwp_ad_g_m2;

% Same ratio with the droplet-mode LWP in the numerator. The wedge is scaled by
% the droplet share of the measured LWP, since the line was fit to the total LWC.
if has_droplet_lwc
    droplet_valid = isfinite(lwc_droplet) & isfinite(z_m);
    lwp_droplet_measured = abs(trapz(z_m(droplet_valid), lwc_droplet(droplet_valid)));   % g/m^2
    droplet_share = lwp_droplet_measured ./ max(lwp_measured, eps);
    adiab.lwp_droplet_g_m2 = lwp_droplet_measured + lwp_wedge .* droplet_share;           % g/m^2
    adiab.f_ad_lwp_droplet = adiab.lwp_droplet_g_m2 ./ adiab.lwp_ad_g_m2;
end

% The lowest sample taken as the base, with no correction of any kind
adiab.lwp_ad_lowest_sample_g_m2 = 0.5 .* Gamma_ad_g_m3_per_m .* h_sampled_m.^2;
adiab.f_ad_lwp_lowest_sample    = lwp_measured ./ adiab.lwp_ad_lowest_sample_g_m2;

end



%% ------------------------------------------------------------------------
%  Thermodynamics. Constants and formulas are copied from
%  lasp-CU-paper-3/plot_per_cloud_rmse_kfold.py -- keep them in sync.
%  ------------------------------------------------------------------------

function con = thermo_constants()

con.G_GRAV  = 9.81;                 % m/s^2
con.R_D     = 287.04;               % J/(kg K)   dry air
con.R_V     = 461.5;                % J/(kg K)   water vapor
con.EPSILON = con.R_D / con.R_V;    % ~0.622
con.C_P     = 1005.0;               % J/(kg K)   dry air at constant pressure
con.L_V     = 2.50e6;               % J/kg       latent heat of vaporization (T ~ 285 K)

end



function e_s = saturation_vapor_pressure_Pa(T_K)
% Bolton (1980): e_s [Pa] as a function of T [K]

T_C = T_K - 273.15;
e_s = 611.2 .* exp(17.67 .* T_C ./ (T_C + 243.5));

end



function q_s = saturation_specific_humidity(T_K, P_Pa)
% Saturation specific humidity [kg/kg]

con = thermo_constants();
e_s = saturation_vapor_pressure_Pa(T_K);
q_s = con.EPSILON .* e_s ./ (P_Pa - (1 - con.EPSILON) .* e_s);

end



function Gamma_m = moist_adiabatic_lapse_rate(T_K, P_Pa)
% Moist adiabatic lapse rate Gamma_m [K/m] at the local state (T, p).
% Randall (2009) eq. 36:
%   Gamma_m = Gamma_d * (1 + L q_s/(R_d T)) / (1 + L^2 q_s/(c_p R_v T^2))
% written with EPSILON/R_D = 1/R_V in the denominator.

con = thermo_constants();
q_s = saturation_specific_humidity(T_K, P_Pa);
num = 1 + con.L_V .* q_s ./ (con.R_D .* T_K);
den = con.C_P + (con.L_V.^2) .* q_s .* con.EPSILON ./ (con.R_D .* T_K.^2);
Gamma_m = con.G_GRAV .* num ./ den;

end



function Gamma_ad = gamma_ad_kg_per_m4(T_K, P_Pa)
% Adiabatic LWC lapse rate Gamma_ad [kg/m^3 per m] at the local state (T, p)
%   Gamma_ad = rho_air * (c_p / L_v) * (Gamma_d - Gamma_m)

con = thermo_constants();
rho_air  = P_Pa ./ (con.R_D .* T_K);           % kg/m^3
Gamma_d  = con.G_GRAV ./ con.C_P;              % K/m
Gamma_m  = moist_adiabatic_lapse_rate(T_K, P_Pa);
Gamma_ad = rho_air .* con.C_P ./ con.L_V .* (Gamma_d - Gamma_m);

end



function adiab = nan_output(z_m, options)
% Every field populated with NaN so callers can always index the structure

n = numel(z_m);

adiab.T_mean_C                   = NaN;
adiab.P_mean_hPa                 = NaN;
adiab.Gamma_m_K_per_km           = NaN;
adiab.Gamma_ad_kg_per_m4         = NaN;
adiab.Gamma_ad_g_m3_per_km       = NaN;
adiab.base_method                = options.base_method;
adiab.base_source                = '';
adiab.n_fit_samples              = NaN;
adiab.fit_slope_g_m3_per_km      = NaN;
adiab.fit_slope_over_gamma_ad    = NaN;
adiab.lwc_fit_at_lowest_g_m3     = NaN;
adiab.z_lowest_sample_m          = NaN;
adiab.dz_base_m                  = NaN;
adiab.z_cb_m                     = NaN;
adiab.z_top_m                    = NaN;
adiab.h_m                        = NaN;
adiab.z_star                     = nan(1, n);
adiab.q_ad_g_m3                  = nan(1, n);
adiab.f_ad_z                     = nan(1, n);
adiab.f_ad_profile               = NaN;
adiab.f_ad_z_droplet             = nan(1, n);
adiab.f_ad_profile_droplet       = NaN;
adiab.lwp_measured_g_m2          = NaN;
adiab.lwp_wedge_g_m2             = NaN;
adiab.lwp_obs_g_m2               = NaN;
adiab.lwp_ad_g_m2                = NaN;
adiab.f_ad_lwp                   = NaN;
adiab.lwp_droplet_g_m2           = NaN;
adiab.f_ad_lwp_droplet           = NaN;
adiab.h_lowest_sample_m          = NaN;
adiab.lwp_ad_lowest_sample_g_m2  = NaN;
adiab.f_ad_lwp_lowest_sample     = NaN;
adiab.lwc_threshold_g_m3         = options.lwc_threshold_g_m3;
adiab.fit_fraction               = options.fit_fraction;
adiab.slope_limits_over_gamma_ad = options.slope_limits_over_gamma_ad;
adiab.zstar_min                  = options.zstar_min;

end
