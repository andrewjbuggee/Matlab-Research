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
%   are averaged over the cloud depth (trapezoid weighting in altitude) and
%   Gamma_ad is evaluated once at that mean state.
%
%   CLOUD BASE. Both methods depend on z_cb, and f_ad_lwp scales as 1/h^2.
%   The profiles handed to this function are truncated where LWC and Nc
%   exceed the profile finder's thresholds, so the lowest sample is NOT the
%   cloud base: its LWC is already 0.05-0.13 g/m^3, which at Gamma_ad ~ 2
%   g/m^3/km puts the true base 25-60 m lower. Using the lowest sample as
%   z_cb inflates f_ad by (1 + dz/h)^2, i.e. 15-35% for 250-400 m clouds,
%   and makes near-adiabatic clouds appear super-adiabatic. By default the
%   base is therefore EXTRAPOLATED: a line is fit to LWC(z) over the lowest
%   third of the profile, its value at the lowest sample LWC_0 is taken, and
%       z_cb = z_lowest - LWC_0 / Gamma_ad
%   i.e. the unsampled layer between the true base and the first sample is
%   assumed adiabatic. This is the smallest correction consistent with the
%   data (a sub-adiabatic base layer would put the base lower still), so any
%   residual bias is still toward high f_ad. The adiabatic wedge below the
%   first sample, (1/2) LWC_0 dz, is added to LWP_obs for consistency.
%   The uncorrected values are kept as *_lowest_sample fields; they are what
%   the Python synthetic-data pipeline currently computes, whose profiles
%   inherit the same truncated base from the in-situ data.
%
%   Wood does not give a formula for Gamma_ad. It follows from conservation
%   of moist static energy along a saturated adiabat,
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
%   same expression appears in Albrecht et al. (1990, GRL 17, 89) and
%   Brenguier (1991, JAS 48, 264).
%
%   The formulas AND the constants are identical to the Python
%   implementation in
%       lasp-CU-paper-3/plot_per_cloud_rmse_kfold.py
%   (saturation_vapor_pressure_Pa, saturation_specific_humidity,
%   moist_adiabatic_lapse_rate, gamma_ad_kg_per_m4, wood_adiabaticity) so
%   the in-situ and synthetic pipelines are directly comparable. If one is
%   changed, change the other.
%
% INPUTS:
%   altitude_m - altitude of each sample [m], in sampling order (a profile
%                may be flown climbing or descending; either is fine)
%   lwc_g_m3   - liquid water content of each sample [g/m^3]
%   temp_C     - ambient air temperature of each sample [deg C]
%   pres_hPa   - ambient static pressure of each sample [hPa]
%   zstar_min  - (optional, default 0.1) lower bound on z* used when
%                averaging f_ad(z) into the scalar f_ad_profile
%   base_method - (optional, default 'extrapolate') how z_cb is defined:
%                'extrapolate'   z_cb = z_lowest - LWC_0/Gamma_ad (see above)
%                'lowest_sample' z_cb = lowest sample altitude (Python pipeline)
%
% OUTPUT:
%   adiab - structure with fields
%       .T_mean_C              mean in-cloud temperature used for Gamma_ad [deg C]
%       .P_mean_hPa            mean in-cloud pressure used for Gamma_ad [hPa]
%       .Gamma_m_K_per_km      moist adiabatic lapse rate at (T_mean, P_mean) [K/km]
%       .Gamma_ad_kg_per_m4    adiabatic LWC lapse rate [kg/m^3 per m]
%       .Gamma_ad_g_m3_per_km  same, in the units people quote [g/m^3 per km]
%       .base_method           'extrapolate' or 'lowest_sample'
%       .z_lowest_sample_m     altitude of the lowest sample [m]
%       .lwc_0_g_m3            fitted LWC at the lowest sample [g/m^3]
%       .dz_base_m             z_lowest_sample - z_cb (0 for 'lowest_sample') [m]
%       .z_cb_m                cloud base used [m]
%       .z_top_m               cloud top = highest sample altitude [m]
%       .h_m                   cloud thickness z_top - z_cb [m]
%       .z_star                (z - z_cb)/h for each sample, sampling order [-]
%       .q_ad_g_m3             adiabatic LWC at each sample, sampling order [g/m^3]
%       .f_ad_z                LWC(z)/q_ad(z) at each sample; NaN where q_ad = 0 [-]
%       .f_ad_profile          method (1) scalar: mean of f_ad_z over z* >= zstar_min [-]
%       .lwp_measured_g_m2     |trapz(altitude, LWC)|, same integral as .lwp [g/m^2]
%       .lwp_obs_g_m2          lwp_measured + adiabatic wedge below the first sample [g/m^2]
%       .lwp_ad_g_m2           (1/2) Gamma_ad h^2 [g/m^2]
%       .f_ad_lwp              method (2) scalar: lwp_obs / lwp_ad [-]
%       .h_lowest_sample_m     z_top - z_lowest_sample [m]
%       .lwp_ad_lowest_sample_g_m2   (1/2) Gamma_ad h_lowest_sample^2 [g/m^2]
%       .f_ad_lwp_lowest_sample      lwp_measured / lwp_ad_lowest_sample: the
%                              uncorrected value, identical in definition to
%                              wood_adiabaticity() in the Python pipeline [-]
%       .zstar_min             value used for f_ad_profile [-]
%
%   Every field is NaN when the profile has no vertical extent, has fewer
%   than two valid (T, p) samples, or yields Gamma_ad <= 0.
%
% SEE ALSO:
%   find_verticalProfiles_VOCALS_REx_ver2, find_verticalProfiles_ORACLES,
%   plot_ensemble_adiabaticity
%
% By Andrew John Buggee
%%

function adiab = compute_adiabaticity_Wood2005(altitude_m, lwc_g_m3, temp_C, pres_hPa, zstar_min, base_method)

if nargin < 5 || isempty(zstar_min)
    zstar_min = 0.1;
end

if nargin < 6 || isempty(base_method)
    base_method = 'extrapolate';
end

if ~any(strcmp(base_method, {'extrapolate', 'lowest_sample'}))
    error([newline, 'base_method must be ''extrapolate'' or ''lowest_sample''.', newline])
end

% Work in double row vectors regardless of how the profile was stored
z_m   = double(altitude_m(:))';          % m
lwc   = double(lwc_g_m3(:))';            % g/m^3
T_K   = double(temp_C(:))' + 273.15;     % K
P_Pa  = double(pres_hPa(:))' .* 100;     % Pa

if numel(z_m) ~= numel(lwc) || numel(z_m) ~= numel(T_K) || numel(z_m) ~= numel(P_Pa)
    error([newline, 'altitude, lwc, temp and pres must have the same number of samples.', newline])
end

adiab = nan_output(z_m, zstar_min, base_method);

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
% Cloud base. Fit LWC(z) over the lowest third of the sampled layer (at
% least 3 samples), evaluate the fit at the lowest sample, and extrapolate
% down the adiabat to LWC = 0. The fit smooths the 1-Hz noise of the
% individual samples nearest the base.
% ------------------------------------------------------------------------
lwp_valid = isfinite(lwc) & isfinite(z_m);

z_star_sampled = (z_m - z_lowest_m) ./ h_sampled_m;
in_lower = lwp_valid & z_star_sampled <= 1/3;

if sum(in_lower) < 3
    [~, order] = sort(z_m);
    order = order(lwp_valid(order));
    in_lower = false(size(z_m));
    in_lower(order(1:min(3, numel(order)))) = true;
end

if sum(in_lower) >= 2 && std(z_m(in_lower)) > 0
    p_fit  = polyfit(z_m(in_lower), lwc(in_lower), 1);
    lwc_0  = polyval(p_fit, z_lowest_m);                % g/m^3, fitted LWC at lowest sample
else
    lwc_0  = mean(lwc(in_lower));
end
lwc_0 = max(lwc_0, 0);

adiab.lwc_0_g_m3 = lwc_0;

if strcmp(base_method, 'extrapolate')
    dz_base_m = lwc_0 ./ Gamma_ad_g_m3_per_m;           % m below the lowest sample
else
    dz_base_m = 0;
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
f_ad_z(~(q_ad > 0)) = NaN;                         % cloud base: 0/0

adiab.q_ad_g_m3 = q_ad;
adiab.f_ad_z    = f_ad_z;

use_for_mean = adiab.z_star >= zstar_min & isfinite(f_ad_z);
if any(use_for_mean)
    adiab.f_ad_profile = mean(f_ad_z(use_for_mean));
end

% ------------------------------------------------------------------------
% Method (2): f_ad = LWP_obs / LWP_ad
% lwp_measured is the same integral as the profile's .lwp field (trapz over
% the sampled altitudes, sign removed). The adiabatic wedge between the
% extrapolated base and the first sample, (1/2) lwc_0 dz_base, is added to
% give lwp_obs; it is zero for base_method = 'lowest_sample'.
% ------------------------------------------------------------------------
lwp_measured = abs(trapz(z_m(lwp_valid), lwc(lwp_valid)));                 % g/m^2
lwp_wedge    = 0.5 .* lwc_0 .* dz_base_m;                                   % g/m^2

adiab.lwp_measured_g_m2 = lwp_measured;
adiab.lwp_obs_g_m2      = lwp_measured + lwp_wedge;
adiab.lwp_ad_g_m2       = 0.5 .* Gamma_ad_g_m3_per_m .* h_m.^2;             % g/m^2
adiab.f_ad_lwp          = adiab.lwp_obs_g_m2 ./ adiab.lwp_ad_g_m2;

% Uncorrected variant: lowest sample taken as the base. Same definition as
% wood_adiabaticity() in lasp-CU-paper-3/plot_per_cloud_rmse_kfold.py.
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



function adiab = nan_output(z_m, zstar_min, base_method)
% Every field populated with NaN so callers can always index the structure

n = numel(z_m);

adiab.T_mean_C                   = NaN;
adiab.P_mean_hPa                 = NaN;
adiab.Gamma_m_K_per_km           = NaN;
adiab.Gamma_ad_kg_per_m4         = NaN;
adiab.Gamma_ad_g_m3_per_km       = NaN;
adiab.base_method                = base_method;
adiab.z_lowest_sample_m          = NaN;
adiab.lwc_0_g_m3                 = NaN;
adiab.dz_base_m                  = NaN;
adiab.z_cb_m                     = NaN;
adiab.z_top_m                    = NaN;
adiab.h_m                        = NaN;
adiab.z_star                     = nan(1, n);
adiab.q_ad_g_m3                  = nan(1, n);
adiab.f_ad_z                     = nan(1, n);
adiab.f_ad_profile               = NaN;
adiab.lwp_measured_g_m2          = NaN;
adiab.lwp_obs_g_m2               = NaN;
adiab.lwp_ad_g_m2                = NaN;
adiab.f_ad_lwp                   = NaN;
adiab.h_lowest_sample_m          = NaN;
adiab.lwp_ad_lowest_sample_g_m2  = NaN;
adiab.f_ad_lwp_lowest_sample     = NaN;
adiab.zstar_min                  = zstar_min;

end
