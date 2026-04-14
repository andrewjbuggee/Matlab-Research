%% ORACLES in-situ vertical profile analysis
%
% Analyzes the ensemble of vertical cloud profiles from ORACLES P-3
% in-situ measurements. Loads the saved mat file produced by
% ensemble_vertical_statistics_ORACLES.m and performs microphysical
% analyses analogous to vocalsRex_inSitu_analysis.m.
%
% Analyses performed:
%   1. Separate profiles: drizzling (D > 50 µm present) vs non-drizzling
%   2. Adiabaticity (Wood 2005) for each profile
%   3. Cloud-Top Entrainment Instability (CTEI) parameter
%   4. Turbulent mixing regime (homogeneous vs inhomogeneous):
%        a. LWC vs phase relaxation time (log-log slope)
%        b. Nc / r_v ratio analysis
%   5. Ensemble MEDIAN profiles of r_e, LWC, N_c for non-precipitating
%      clouds vs normalized optical depth
%
% By Andrew John Buggee

clear variables

%% -------------------------------------------------------------------------
%  Physical constants
%  -------------------------------------------------------------------------

g     = 9.81;        % m/s²  - gravitational acceleration
c_pd  = 1005;        % J/(kg·K) - specific heat of dry air at const. pressure
L_v   = 2.501e6;     % J/kg  - latent heat of vaporization (0°C)
R_d   = 287.04;      % J/(kg·K) - specific gas constant, dry air
R_v   = 461.5;       % J/(kg·K) - specific gas constant, water vapor
epsilon = R_d / R_v; % ≈ 0.622 - ratio of molecular weights
rho_L = 1000;        % kg/m³ - density of liquid water
D_v0  = 2.21e-5;     % m²/s  - diffusivity of vapor in air at 0°C, 1013 mb


%% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


%% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

% Search for ensemble profiles mat file (use the most recently created one)
folder_contents = dir(foldername_data);
mat_candidates = {};

for nn = 1:length(folder_contents)
    fname = folder_contents(nn).name;
    if length(fname) > 17 && strcmp(fname(1:17), 'ensemble_profiles') && ...
            strcmp(fname(end-3:end), '.mat')
        mat_candidates{end+1} = fname; %#ok<SAGROW>
    end
end

if isempty(mat_candidates)
    error('No ensemble_profiles mat file found in:\n  %s', foldername_data)
end

% Use the last file listed (typically most recent by name/date)
mat_filename = mat_candidates{end};
disp(['Loading: ', mat_filename])
load([foldername_data, mat_filename], 'ensemble_profiles', 'inputs')

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])


%% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])

% --- Pie / bar summary plot ---
plt_clr_1 = mySavedColors(62, 'fixed');
plt_clr_2 = mySavedColors(62, 'fixed');
figure;
num_driz_nonDriz = [length(idx_no_drizzle), length(idx_drizzle)];
bar(num_driz_nonDriz, 0.5, ...
    'FaceColor', plt_clr_1)
text(1:length(num_driz_nonDriz), num_driz_nonDriz, num2str(num_driz_nonDriz'),...
    'vert','bottom','horiz','center');
set(gca, 'XTickLabel', {'Non-drizzling', 'Drizzling'}, 'Fontsize', 20)
ylabel('Number of profiles', 'fontsize', 20)
title(['ORACLES vertical profiles', newline', 'drizzle LWP threshold ($LWP_{2DS} > $', ...
    num2str(drizzle_LWP_threshold), ' $g/m^{2}$)'], 'interpreter', 'latex',...
    'fontsize', 20)
grid on; grid minor
set(gcf, 'Position', [0 0 600 600])


%% -------------------------------------------------------------------------
%  2. Adiabaticity (Wood 2005) for each profile
%
%  Adiabatic LWC lapse rate (Brenguier et al. 2000):
%    Γ_L = ρ_a * (DALR - SALR)     [g m⁻⁴]
%
%  where DALR = g/c_pd and SALR = g*(1 + Lv*rs/(Rd*T))/(c_pd + Lv²*rs/(Rv*T²))
%
%  Adiabatic LWP: LWP_ad = 0.5 * Γ_L * H²
%  Adiabaticity:  A_d = LWP_meas / LWP_ad
%  -------------------------------------------------------------------------

adiabaticity      = NaN(1, N_profiles);
LWP_measured      = NaN(1, N_profiles);
LWP_adiabatic     = NaN(1, N_profiles);
Gamma_L_profile   = NaN(1, N_profiles);   % adiabatic LWC lapse rate g/m^4

for nn = 1:N_profiles

    [re_p, lwc_p, Nc_p, alt_p, T_K_p, P_Pa_p, rv_p] = orientProfile(ensemble_profiles{nn});

    % Cloud base: index 1 (after orientation to cloud-base-first)
    T_base = T_K_p(1);       % K
    P_base = P_Pa_p(1);      % Pa

    % Saturation vapor pressure (Tetens formula) [Pa]
    e_s_base = 611.2 * exp(17.67 * (T_base - 273.15) / (T_base - 29.65));

    % Saturation mixing ratio at cloud base [kg/kg]
    r_s_base = epsilon * e_s_base / max(P_base - e_s_base, 1);

    % Air density at cloud base [kg/m³]
    rho_a_base = P_base / (R_d * T_base);

    % Dry adiabatic lapse rate [K/m]
    DALR = g / c_pd;

    % Moist/saturated adiabatic lapse rate at cloud base [K/m]
    SALR = g * (1 + L_v * r_s_base / (R_d * T_base)) / ...
        (c_pd + L_v^2 * r_s_base / (R_v * T_base^2));

    % Adiabatic LWC lapse rate [kg/m^4], converted to [g/m^4]
    Gamma_L = rho_a_base * (DALR - SALR) * 1000;    % g/m^4
    Gamma_L_profile(nn) = Gamma_L;

    % Height above cloud base [m]
    H_total = abs(alt_p(end) - alt_p(1));           % cloud geometric depth [m]
    z_cb    = alt_p - alt_p(1);                     % height above cloud base

    % Adiabatic LWC profile [g/m³]
    LWC_ad_profile = max(Gamma_L * z_cb, 0);

    % Measured and adiabatic LWP [g/m²]
    LWP_measured(nn)  = trapz(alt_p, lwc_p);        % g/m²  (always positive)
    LWP_adiabatic(nn) = 0.5 * Gamma_L * H_total^2; % g/m²

    if LWP_adiabatic(nn) > 0
        adiabaticity(nn) = LWP_measured(nn) / LWP_adiabatic(nn);
    end

end

% --- Plot histogram of adiabaticity for drizzle vs non-drizzle ---
figure1 = figure;
hold on
histogram(adiabaticity(~is_drizzle), 'BinWidth', 0.05, ...
    'FaceColor', plt_clr_1, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'Normalization', 'probability')
histogram(adiabaticity(is_drizzle), 'BinWidth', 0.05, ...
    'FaceColor', plt_clr_2, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'Normalization', 'probability')
xline(1, 'k--', 'LineWidth', 1.5)
xlabel('Adiabaticity  $A_d = LWP / LWP_{ad}$', 'Interpreter', 'latex')
ylabel('Relative frequency')
title('Cloud Adiabaticity — ORACLES', 'Interpreter', 'latex')
legend({'Non-drizzling', 'Drizzling', 'Adiabatic'}, 'Location', 'best')
grid on; grid minor
set(gcf, 'Position', [0 0 650 500])


% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'Histogram of Adiabaticity from ORACLES in-situ measurements.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs, 'Histogram of Adiabaticity from ORACLES in-situ measurements.png'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------



% --- Adiabaticity vs LWP scatter ---
figure; hold on
scatter(LWP_measured(~is_drizzle), adiabaticity(~is_drizzle), 30, ...
    plt_clr_1, 'filled', 'MarkerFaceAlpha', 0.5)
scatter(LWP_measured(is_drizzle), adiabaticity(is_drizzle), 30, ...
    plt_clr_2, 'filled', 'MarkerFaceAlpha', 0.5)
yline(1, 'k--', 'LineWidth', 1.5)
xlabel('Measured LWP (g/m²)', 'Interpreter', 'latex')
ylabel('Adiabaticity  $A_d$', 'Interpreter', 'latex')
title('Adiabaticity vs LWP — ORACLES', 'Interpreter', 'latex')
legend({'Non-drizzling', 'Drizzling', 'Adiabatic ($A_d=1$)'}, ...
    'Location', 'best', 'Interpreter', 'latex')
grid on; grid minor
set(gcf, 'Position', [0 0 650 500])






%% -------------------------------------------------------------------------
%  3. Cloud-Top Entrainment Instability (CTEI)
%
%  CTEI criterion (Randall 1980, Deardorff 1980, Kuo & Schubert 1988):
%    κ = c_p * ΔT / (L_v * |Δq_s|)    cloud is buoyancy unstable when κ > κ_c
%
%  where ΔT = T_FT - T_CT > 0 (warmer free troposphere above inversion)
%        Δq_s = q_s,FT - q_s,CT < 0 (drier free troposphere)
%        κ_c ≈ 0.23  (critical value, Kuo & Schubert 1988)
%
%  Since in-situ profiles do not include above-cloud measurements, we
%  estimate the inversion jump (ΔT) from the departure of the cloud-top
%  temperature from the moist adiabat predicted from cloud base. The
%  excess warming above the moist adiabat at cloud top is taken as a
%  lower-bound estimate of the inversion strength.
%
%  Liquid-water potential temperature (θ_L) is also computed as a
%  thermodynamic tracer; large θ_L at cloud top indicates strong inversion.
%  -------------------------------------------------------------------------

kappa_CTEI  = NaN(1, N_profiles);  % CTEI instability parameter κ
theta_L_top = NaN(1, N_profiles);  % liquid water potential temperature at CT
theta_v_top = NaN(1, N_profiles);  % virtual potential temperature at CT
DeltaT_inv  = NaN(1, N_profiles);  % estimated inversion strength [K]

kappa_c = 0.23;   % critical value for CTEI (Kuo & Schubert 1988)

P0 = 1e5;   % reference pressure for potential temperature [Pa]

for nn = 1:N_profiles

    [~, lwc_p, ~, alt_p, T_K_p, P_Pa_p, ~] = orientProfile(ensemble_profiles{nn});

    % Cloud-base and cloud-top indices (after orientation: base = 1, top = end)
    T_base = T_K_p(1);
    P_base = P_Pa_p(1);
    T_top  = T_K_p(end);
    P_top  = P_Pa_p(end);

    % --- Liquid-water potential temperature at cloud top ---
    % θ_L ≈ θ * (1 - L_v*q_L / (c_pd*T))
    % where θ = T*(P0/P)^(Rd/cpd) and q_L = LWC/ρ_a [kg/kg]
    theta_top = T_top * (P0 / P_top)^(R_d / c_pd);
    rho_a_top = P_top / (R_d * T_top);
    q_L_top   = lwc_p(end) / (rho_a_top * 1000);   % kg/kg (LWC in g/m³)
    theta_L_top(nn) = theta_top * (1 - L_v * q_L_top / (c_pd * T_top));

    % --- Virtual potential temperature at cloud top (assume saturated) ---
    e_s_top = 611.2 * exp(17.67 * (T_top - 273.15) / (T_top - 29.65));
    q_s_top = epsilon * e_s_top / max(P_top - e_s_top, 1);   % kg/kg
    q_v_top = q_s_top;   % saturated inside cloud
    theta_v_top(nn) = theta_top * (1 + 0.61 * q_v_top - q_L_top);

    % --- Estimate inversion strength from moist-adiabat departure ---
    % Predict T at cloud top if lifting from cloud base were purely moist adiabatic
    e_s_base  = 611.2 * exp(17.67 * (T_base - 273.15) / (T_base - 29.65));
    r_s_base  = epsilon * e_s_base / max(P_base - e_s_base, 1);
    rho_a_base = P_base / (R_d * T_base);
    DALR = g / c_pd;
    SALR = g * (1 + L_v * r_s_base / (R_d * T_base)) / ...
        (c_pd + L_v^2 * r_s_base / (R_v * T_base^2));
    H_total = abs(alt_p(end) - alt_p(1));   % cloud depth [m]
    T_moist_adiabat_top = T_base - SALR * H_total;   % predicted T at cloud top

    % Departure: measured T_top is warmer than moist adiabat → inversion entraining
    DeltaT_inv(nn) = T_top - T_moist_adiabat_top;   % K; positive = warmer than expected

    % --- CTEI parameter κ ---
    % Use estimated ΔT = DeltaT_inv as a lower-bound inversion jump
    % Δq_s = qs(T_FT) - qs(T_top) where T_FT = T_top + DeltaT_inv
    if DeltaT_inv(nn) > 0
        T_FT_est  = T_top + DeltaT_inv(nn);
        e_s_FT    = 611.2 * exp(17.67 * (T_FT_est - 273.15) / (T_FT_est - 29.65));
        q_s_FT    = epsilon * e_s_FT / max(P_top - e_s_FT, 1);
        Delta_q_s = q_s_top - q_s_FT;   % < 0: FT is drier
        if abs(Delta_q_s) > 1e-6
            kappa_CTEI(nn) = c_pd * DeltaT_inv(nn) / (L_v * abs(Delta_q_s));
        end
    end

end

% --- CTEI plots ---
figure;
subplot(1,2,1)
hold on
histogram(kappa_CTEI(~is_drizzle), 'BinWidth', 0.05, ...
    'FaceColor', plt_clr_1, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'Normalization', 'probability')
histogram(kappa_CTEI(is_drizzle), 'BinWidth', 0.05, ...
    'FaceColor', plt_clr_2, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'Normalization', 'probability')
xline(kappa_c, 'k--', 'LineWidth', 1.5)
xlabel('$\kappa_{CTEI}$', 'Interpreter', 'latex')
ylabel('Relative frequency')
title('CTEI parameter $\kappa$', 'Interpreter', 'latex')
legend({'Non-drizzling', 'Drizzling', ['$\kappa_c = ', num2str(kappa_c), '$']}, ...
    'Location', 'best', 'Interpreter', 'latex')
grid on; grid minor

subplot(1,2,2)
hold on
scatter(LWP_measured(~is_drizzle), kappa_CTEI(~is_drizzle), 30, ...
    plt_clr_1, 'filled', 'MarkerFaceAlpha', 0.5)
scatter(LWP_measured(is_drizzle), kappa_CTEI(is_drizzle), 30, ...
    plt_clr_2, 'filled', 'MarkerFaceAlpha', 0.5)
yline(kappa_c, 'k--', 'LineWidth', 1.5)
xlabel('Measured LWP (g/m²)', 'Interpreter', 'latex')
ylabel('$\kappa_{CTEI}$', 'Interpreter', 'latex')
title('CTEI vs LWP', 'Interpreter', 'latex')
legend({'Non-drizzling', 'Drizzling'}, 'Location', 'best')
grid on; grid minor

set(gcf, 'Position', [0 0 1100 500])


%% -------------------------------------------------------------------------
%  4a. Turbulent mixing regime: LWC vs phase relaxation time (τ_phase)
%
%  Phase relaxation time (Korolev & Isaac 2000, Lehmann et al. 2009):
%    τ_phase = ρ_L / (4π * D_v * N_c * r̄)
%
%  where D_v = D_v0 * (T/273.15)^1.94 * (1013 mb / P_mb) [m²/s]
%        N_c in m⁻³, r̄ = mean radius ≈ r_v [m]
%
%  On a log-log plot, slope of LWC vs τ_phase:
%    homogeneous mixing   → slope ≈ -2/3 to -1
%    extreme inhomogeneous → slope ≈ 0
%  -------------------------------------------------------------------------

% Collect all in-cloud LWC and τ_phase values
all_LWC_nodrizzle    = [];
all_tau_phase_nodrizzle = [];
all_LWC_drizzle      = [];
all_tau_phase_drizzle   = [];

for nn = 1:N_profiles

    [~, lwc_p, Nc_p, ~, T_K_p, P_Pa_p, rv_p] = orientProfile(ensemble_profiles{nn});

    % Diffusivity of vapor [m²/s] at each level
    P_mb_p = P_Pa_p / 100;
    D_v_p  = D_v0 * (T_K_p / 273.15).^1.94 .* (1013 ./ P_mb_p);

    % Number concentration in m⁻³ (convert from cm⁻³)
    Nc_m3 = Nc_p * 1e6;   % m⁻³

    % Mean radius [m]; use r_v (mean volume radius) in µm → m
    r_mean_m = rv_p * 1e-6;   % µm → m

    % Phase relaxation time [s]
    denom = 4 * pi * D_v_p .* Nc_m3 .* r_mean_m;
    valid = denom > 0 & lwc_p > 0 & Nc_p > 0 & rv_p > 0;
    tau_phase_p = NaN(size(lwc_p));
    tau_phase_p(valid) = rho_L ./ denom(valid);

    % Accumulate
    if is_drizzle(nn)
        all_LWC_drizzle    = [all_LWC_drizzle,    lwc_p(valid)];
        all_tau_phase_drizzle = [all_tau_phase_drizzle, tau_phase_p(valid)]; %#ok<*AGROW>
    else
        all_LWC_nodrizzle    = [all_LWC_nodrizzle,    lwc_p(valid)];
        all_tau_phase_nodrizzle = [all_tau_phase_nodrizzle, tau_phase_p(valid)];
    end

end

% --- Estimate mixing regime slope via log-log linear regression ---
coeff_nd = polyfit(log10(all_tau_phase_nodrizzle), log10(all_LWC_nodrizzle), 1);
coeff_d  = polyfit(log10(all_tau_phase_drizzle),   log10(all_LWC_drizzle),   1);

tau_fit  = logspace(log10(min([all_tau_phase_nodrizzle, all_tau_phase_drizzle])), ...
    log10(max([all_tau_phase_nodrizzle, all_tau_phase_drizzle])), 100);
LWC_fit_nd = 10.^polyval(coeff_nd, log10(tau_fit));
LWC_fit_d  = 10.^polyval(coeff_d,  log10(tau_fit));

% --- Plot LWC vs τ_phase ---
figure;
subplot(1,2,1)
hold on
scatter(all_tau_phase_nodrizzle, all_LWC_nodrizzle, 2, ...
    plt_clr_1, 'filled', 'MarkerFaceAlpha', 0.05)
plot(tau_fit, LWC_fit_nd, '-', 'Color', plt_clr_1, 'LineWidth', 2)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on; grid minor
xlabel('Phase relaxation time $\tau_{ph}$ (s)', 'Interpreter', 'latex')
ylabel('LWC (g/m³)', 'Interpreter', 'latex')
title(['Non-drizzling  |  slope = ', sprintf('%.2f', coeff_nd(1))], 'Interpreter', 'latex')

subplot(1,2,2)
hold on
scatter(all_tau_phase_drizzle, all_LWC_drizzle, 2, ...
    plt_clr_2, 'filled', 'MarkerFaceAlpha', 0.05)
plot(tau_fit, LWC_fit_d, '-', 'Color', plt_clr_2, 'LineWidth', 2)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on; grid minor
xlabel('Phase relaxation time $\tau_{ph}$ (s)', 'Interpreter', 'latex')
ylabel('LWC (g/m³)', 'Interpreter', 'latex')
title(['Drizzling  |  slope = ', sprintf('%.2f', coeff_d(1))], 'Interpreter', 'latex')

sgtitle({'LWC vs Phase Relaxation Time — ORACLES', ...
    'Mixing regime: homogeneous $\approx -1$, inhomogeneous $\approx 0$'}, ...
    'Interpreter', 'latex')
set(gcf, 'Position', [0 0 1100 500])


%% -------------------------------------------------------------------------
%  4b. Turbulent mixing regime: N_c / r_v ratio
%
%  Under homogeneous mixing: both N_c and r_v change with dilution, but
%    N_c / r_v^3 ∝ LWC (both decrease together → Nc/rv changes)
%  Under inhomogeneous mixing: r_v is conserved, only N_c decreases,
%    so N_c / r_v remains proportional to N_c alone
%
%  A plot of N_c vs r_v on log-log axes with reference lines for
%  constant LWC distinguishes the two regimes (Gerber 2000,
%  Pawlowska et al. 2006):
%    homogeneous:     data scatter along constant-LWC lines
%    inhomogeneous:   data align vertically (constant r_v, varying N_c)
%  -------------------------------------------------------------------------

% Collect N_c and r_v across all points
all_Nc_nd = [];  all_rv_nd = [];  all_lwc_nd_b = [];
all_Nc_d  = [];  all_rv_d  = [];  all_lwc_d_b  = [];

for nn = 1:N_profiles

    [~, lwc_p, Nc_p, ~, ~, ~, rv_p] = orientProfile(ensemble_profiles{nn});
    valid = Nc_p > 0 & rv_p > 0 & lwc_p > 0;

    if is_drizzle(nn)
        all_Nc_d     = [all_Nc_d,     Nc_p(valid)];
        all_rv_d     = [all_rv_d,     rv_p(valid)];
        all_lwc_d_b  = [all_lwc_d_b,  lwc_p(valid)];
    else
        all_Nc_nd    = [all_Nc_nd,    Nc_p(valid)];
        all_rv_nd    = [all_rv_nd,    rv_p(valid)];
        all_lwc_nd_b = [all_lwc_nd_b, lwc_p(valid)];
    end

end

% Reference lines of constant LWC: LWC = (4π/3)*ρ_L*N_c*r_v^3
%   → N_c = LWC / [(4π/3)*ρ_L*r_v^3]   (LWC in g/m³, r_v in µm → m)
rv_line  = logspace(0, 2.5, 200);  % µm
LWC_lines = [0.05, 0.1, 0.2, 0.4, 0.6];   % g/m³

figure;
subplot(1,2,1)
hold on
scatter(all_rv_nd, all_Nc_nd, 2, all_lwc_nd_b, 'filled', 'MarkerFaceAlpha', 0.1)
colormap(gca, 'turbo')
cb = colorbar; cb.Label.String = 'LWC (g/m³)';
for ll = 1:length(LWC_lines)
    Nc_ref = LWC_lines(ll) ./ ((4*pi/3) * (rho_L*1e-3) * (rv_line*1e-6).^3 * 1e6);
    plot(rv_line, Nc_ref, 'k-', 'LineWidth', 0.8)
    text(rv_line(end), Nc_ref(end), [num2str(LWC_lines(ll)),' g/m³'], ...
        'FontSize', 7, 'HorizontalAlignment', 'right')
end
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$r_v$ (µm)', 'Interpreter', 'latex')
ylabel('$N_c$ (cm$^{-3}$)', 'Interpreter', 'latex')
title('Non-drizzling', 'Interpreter', 'latex')
grid on; grid minor

subplot(1,2,2)
hold on
scatter(all_rv_d, all_Nc_d, 2, all_lwc_d_b, 'filled', 'MarkerFaceAlpha', 0.1)
colormap(gca, 'turbo')
cb = colorbar; cb.Label.String = 'LWC (g/m³)';
for ll = 1:length(LWC_lines)
    Nc_ref = LWC_lines(ll) ./ ((4*pi/3) * (rho_L*1e-3) * (rv_line*1e-6).^3 * 1e6);
    plot(rv_line, Nc_ref, 'k-', 'LineWidth', 0.8)
    text(rv_line(end), Nc_ref(end), [num2str(LWC_lines(ll)),' g/m³'], ...
        'FontSize', 7, 'HorizontalAlignment', 'right')
end
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$r_v$ (µm)', 'Interpreter', 'latex')
ylabel('$N_c$ (cm$^{-3}$)', 'Interpreter', 'latex')
title('Drizzling', 'Interpreter', 'latex')
grid on; grid minor

sgtitle({'$N_c$ vs $r_v$ — ORACLES  |  Mixing regime diagnostics', ...
    'Reference lines = constant LWC; vertical scatter → inhomogeneous'}, ...
    'Interpreter', 'latex')
set(gcf, 'Position', [0 0 1200 550])



%% -------------------------------------------------------------------------
%  5 & 6. Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd  = zeros(n_bins,1);  re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd = zeros(n_bins,1);  lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd  = zeros(n_bins,1);  Nc_iqr_nd  = zeros(n_bins,1);

re_med_d   = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d  = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d   = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_med_nd(bb)  = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb) = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)  = median(vertSeg_nd{bb,3}, 'omitnan');
    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_med_d(bb)   = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});
end

% --- Plot ensemble median profiles for non-drizzling clouds ---
figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 22;
ttl_fnt = 26;

% r_e
subplot(1,3,1)
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, max(re_med_nd)*1.3])

% LWC
subplot(1,3,2)
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['Median profiles — Non-drizzling ORACLES  |  ', ...
    num2str(length(idx_no_drizzle)), ' profiles'], 'Interpreter', 'latex', 'FontSize', ttl_fnt)

% N_c
subplot(1,3,3)
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

set(gcf, 'Position', [0 0 1200 600])






% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'Position',[0.185000000000001 0.11 0.21340579710145 0.815]);
hold(axes1,'on');

subplot(1,3,1)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.00447158452278369 0.806666666666669 0.14136174881055 0.111660079616977])







subplot(1,3,2)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title('ORACLES median profiles: drizzling vs non-drizzling', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)



subplot(1,3,3)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

set(gcf, 'Position', [0 0 1200 600])




% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% % save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
% end
% saveas(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle',...
%     '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------



%% -------------------------------------------------------------------------
%  5 & 6. Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version).  Take 2
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Create a PDF object at each level in the cloud and fit a distribution to this PDF

% store the refection of each null hypothesis and the p-value for each
% chi-squared test

% re_reject_normal_nd = zeros(1, size(vertSeg_nd,1));
% re_p_normal_nd = zeros(1, size(vertSeg_nd,1));

re_reject_lognormal_nd = zeros(1, size(vertSeg_nd,1));
re_p_lognormal_nd = zeros(1, size(vertSeg_nd,1));

% re_reject_gamma_nd = zeros(1, size(vertSeg_nd,1));
% re_p_gamma_nd = zeros(1, size(vertSeg_nd,1));


re_reject_normal_d = zeros(1, size(vertSeg_d,1));
re_p_normal_d = zeros(1, size(vertSeg_d,1));

re_reject_lognormal_d = zeros(1, size(vertSeg_d,1));
re_p_lognormal_d = zeros(1, size(vertSeg_d,1));

re_reject_gamma_d = zeros(1, size(vertSeg_d,1));
re_p_gamma_d = zeros(1, size(vertSeg_d,1));






lwc_reject_normal_nd = zeros(1, size(vertSeg_nd,1));
lwc_p_normal_nd = zeros(1, size(vertSeg_nd,1));

% lwc_reject_lognormal = zeros(1, size(vertSeg_nd,1));
% lwc_p_lognormal = zeros(1, size(vertSeg_nd,1));
%
% lwc_reject_gamma = zeros(1, size(vertSeg_nd,1));
% lwc_p_gamma = zeros(1, size(vertSeg_nd,1));


lwc_reject_normal_d = zeros(1, size(vertSeg_d,1));
lwc_p_normal_d = zeros(1, size(vertSeg_d,1));

lwc_reject_lognormal_d = zeros(1, size(vertSeg_d,1));
lwc_p_lognormal_d = zeros(1, size(vertSeg_d,1));

lwc_reject_gamma_d = zeros(1, size(vertSeg_d,1));
lwc_p_gamma_d = zeros(1, size(vertSeg_d,1));





% Nc_reject_normal = zeros(1, size(vertSeg_nd,1));
% Nc_p_normal = zeros(1, size(vertSeg_nd,1));
%
% Nc_reject_lognormal = zeros(1, size(vertSeg_nd,1));
% Nc_p_lognormal = zeros(1, size(vertSeg_nd,1));

Nc_reject_gamma_nd = zeros(1, size(vertSeg_nd,1));
Nc_p_gamma_nd = zeros(1, size(vertSeg_nd,1));


Nc_reject_normal_d = zeros(1, size(vertSeg_d,1));
Nc_p_normal_d = zeros(1, size(vertSeg_d,1));

Nc_reject_lognormal_d = zeros(1, size(vertSeg_d,1));
Nc_p_lognormal_d = zeros(1, size(vertSeg_d,1));

Nc_reject_gamma_d = zeros(1, size(vertSeg_d,1));
Nc_p_gamma_d = zeros(1, size(vertSeg_d,1));







for bb = 1:size(vertSeg_nd, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    % re_fit_normal(bb) = fitdist(vertSeg_nd{bb,1}', 'normal');
    % [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertSeg_nd{bb,1}', 'CDF',re_fit_normal(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertSeg_nd{bb,1}', 'lognormal');
    [re_reject_lognormal_nd(bb), re_p_lognormal_nd(bb)] = chi2gof(vertSeg_nd{bb,1}', 'CDF',re_fit_lognormal(bb));

    % fit the effective radius data to a gamma distribution
    % re_fit_gamma(bb) = fitdist(vertSeg_nd{bb,1}', 'gamma');
    % [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertSeg_nd{bb,1}', 'CDF', re_fit_gamma(bb));

    
    % fit the effective radius data to a normal distribution
    re_fit_normal_d(bb) = fitdist(vertSeg_d{bb,1}', 'normal');
    [re_reject_normal_d(bb), re_p_normal_d(bb)] = chi2gof(vertSeg_d{bb,1}', 'CDF',re_fit_normal_d(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal_d(bb) = fitdist(vertSeg_d{bb,1}', 'lognormal');
    [re_reject_lognormal_n(bb), re_p_lognormal_d(bb)] = chi2gof(vertSeg_d{bb,1}', 'CDF',re_fit_lognormal_d(bb));

    % fit the effective radius data to a gamma distribution
    re_fit_gamma_d(bb) = fitdist(vertSeg_d{bb,1}', 'gamma');
    [re_reject_gamma_d(bb), re_p_gamma_d(bb)] = chi2gof(vertSeg_d{bb,1}', 'CDF', re_fit_gamma_d(bb));







    % -------------------------------------------
    % ------- LIQUID WATER CONTENT FITTING ------
    % -------------------------------------------


    % fit the liquid water content data to a normal distribution
    lwc_fit_normal(bb) = fitdist(vertSeg_nd{bb,2}', 'normal');
    [lwc_reject_normal_nd(bb), lwc_p_normal_nd(bb)] = chi2gof(vertSeg_nd{bb,2}', 'CDF',lwc_fit_normal(bb));

    % % fit the liquid water content data to a log-normal distribution
    % lwc_fit_lognormal(bb) = fitdist(vertSeg_nd{bb,2}', 'lognormal');
    % [lwc_reject_lognormal(bb), lwc_p_lognormal(bb)] = chi2gof(vertSeg_nd{bb,2}', 'CDF',lwc_fit_lognormal(bb));
    %
    % % fit the liquid water content data to a gamma distribution
    % lwc_fit_gamma(bb) = fitdist(vertSeg_nd{bb,2}', 'gamma');
    % [lwc_reject_gamma(bb), lwc_p_gamma(bb)] = chi2gof(vertSeg_nd{bb,2}', 'CDF', lwc_fit_gamma(bb));

    
    % fit the liquid water content data to a normal distribution
    lwc_fit_normal_d(bb) = fitdist(vertSeg_d{bb,2}', 'normal');
    [lwc_reject_normal_d(bb), lwc_p_normal_d(bb)] = chi2gof(vertSeg_d{bb,2}', 'CDF', lwc_fit_normal_d(bb));

    % fit the liquid water content data to a log-normal distribution
    lwc_fit_lognormal_d(bb) = fitdist(vertSeg_d{bb,2}', 'lognormal');
    [lwc_reject_lognormal_n(bb), lwc_p_lognormal_d(bb)] = chi2gof(vertSeg_d{bb,2}', 'CDF',lwc_fit_lognormal_d(bb));

    % fit the liquid water content data to a gamma distribution
    lwc_fit_gamma_d(bb) = fitdist(vertSeg_d{bb,2}', 'gamma');
    [lwc_reject_gamma_d(bb), lwc_p_gamma_d(bb)] = chi2gof(vertSeg_d{bb,2}', 'CDF', lwc_fit_gamma_d(bb));







    % -------------------------------------------
    % ------- NUMBER CONCENTRATION FITTING ------
    % -------------------------------------------


    % % fit the number concentration data to a normal distribution
    % Nc_fit_normal(bb) = fitdist(vertSeg_nd{bb,3}', 'normal');
    % [Nc_reject_normal(bb), Nc_p_normal(bb)] = chi2gof(vertSeg_nd{bb,3}', 'CDF',Nc_fit_normal(bb));
    %
    % % fit the number concentration content data to a log-normal distribution
    % Nc_fit_lognormal(bb) = fitdist(vertSeg_nd{bb,3}', 'lognormal');
    % [Nc_reject_lognormal(bb), Nc_p_lognormal(bb)] = chi2gof(vertSeg_nd{bb,3}', 'CDF',Nc_fit_lognormal(bb));

    % fit the number concentration content data to a gamma distribution
    Nc_fit_gamma_nd(bb) = fitdist(vertSeg_nd{bb,3}', 'gamma');
    [Nc_reject_gamma_nd(bb), Nc_p_gamma_nd(bb)] = chi2gof(vertSeg_nd{bb,3}', 'CDF', Nc_fit_gamma_nd(bb));


    % fit the number concentration data to a normal distribution
    Nc_fit_normal_d(bb) = fitdist(vertSeg_d{bb,3}', 'normal');
    [Nc_reject_normal_d(bb), Nc_p_normal_d(bb)] = chi2gof(vertSeg_d{bb,3}', 'CDF',Nc_fit_normal_d(bb));

    % fit the number concentration content data to a log-normal distribution
    Nc_fit_lognormal_d(bb) = fitdist(vertSeg_d{bb,3}', 'lognormal');
    [Nc_reject_lognormal_d(bb), Nc_p_lognormal_d(bb)] = chi2gof(vertSeg_d{bb,3}', 'CDF',Nc_fit_lognormal_d(bb));

    % fit the number concentration content data to a gamma distribution
    Nc_fit_gamma_d(bb) = fitdist(vertSeg_d{bb,3}', 'gamma');
    [Nc_reject_gamma_d(bb), Nc_p_gamma_d(bb)] = chi2gof(vertSeg_d{bb,3}', 'CDF', Nc_fit_gamma_d(bb));






end



% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

% bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% % -----------------------------------------------
% % ------- EFFECTIVE DROPLET RADIUS FITTING ------
% % -----------------------------------------------
% [max_re_p, idx_re_p] = max([re_p_normal; re_p_lognormal; re_p_gamma],[], 1);
%
% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]);
% title('r_e best distribution fit'); ylabel('Counts')
%
%
%
% % -------------------------------------------
% % ------- LIQUID WATER CONTENT FITTING ------
% % -------------------------------------------
% [max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1);
%
% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]);
% title('LWC best distribution fit'); ylabel('Counts')
%
%
% % -------------------------------------------
% % ------- NUMBER CONCENTRATION FITTING ------
% % -------------------------------------------
%
% [max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1);
%
% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]);
% title('N_c best distribution fit'); ylabel('Counts')



bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
[max_re_p_d, idx_re_p_d] = max([re_p_normal_d; re_p_lognormal_d; re_p_gamma_d],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p_d==1), sum(idx_re_p_d==2), sum(idx_re_p_d==3)]);
title('r_e best distribution fit (w/ drizzle)'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p_d, idx_lwc_p_d] = max([lwc_p_normal_d; lwc_p_lognormal_d; lwc_p_gamma_d],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p_d==1), sum(idx_lwc_p_d==2), sum(idx_lwc_p_d==3)]);
title('LWC best distribution fit (w/ drizzle)'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p_d, idx_Nc_p_d] = max([Nc_p_normal_d; Nc_p_lognormal_d; Nc_p_gamma_d],[], 1);

figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p_d==1), sum(idx_Nc_p_d==2), sum(idx_Nc_p_d==3)]);
title('N_c best distribution fit (w/ drizzle)'); ylabel('Counts')





% Because of the asymmetry of log-normal distributions, it is indeed common
%  — and arguably more natural — to characterize lognormal distributions using
%  the geometric mean and geometric standard deviation (GSD) rather than
%  the arithmetic equivalents.

% If ln(X) ~ N(μ, σ²), then:
%   Geometric mean (GM): exp(μ) — this equals the median of X, which for a lognormal is a more representative center than the arithmetic mean.
%   Geometric standard deviation (GSD): exp(σ)


% Compute median and IQR for each bin
re_GM_nd  = zeros(n_bins,1);  re_GSD_nd  = zeros(n_bins,1);
lwc_GM_nd = zeros(n_bins,1);  lwc_GSD_nd = zeros(n_bins,1);
Nc_GM_nd  = zeros(n_bins,1);  Nc_GSD_nd  = zeros(n_bins,1);

re_GM_d   = zeros(n_bins,1);  re_GSD_d   = zeros(n_bins,1);
lwc_GM_d  = zeros(n_bins,1);  lwc_GSD_d  = zeros(n_bins,1);
Nc_GM_d   = zeros(n_bins,1);  Nc_GSD_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_GM_nd(bb)  = exp(re_fit_lognormal(bb).mu);         % microns
    lwc_GM_nd(bb) = lwc_fit_lognormal(bb).mu;
    Nc_GM_nd(bb)  = exp(re_fit_lognormal(bb).mu);
    re_GSD_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_GSD_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_GSD_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_GSD_d(bb)   = exp(re_fit_lognormal(bb).mu);         % microns
    lwc_GSD_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_GSD_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_GSD_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_GSD_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_GSD_d(bb)   = iqr(vertSeg_d{bb,3});
end

% --- Plot ensemble median profiles for non-drizzling clouds ---
figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 22;
ttl_fnt = 26;

% r_e
subplot(1,3,1)
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, max(re_med_nd)*1.3])

% LWC
subplot(1,3,2)
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['Median profiles — Non-drizzling ORACLES  |  ', ...
    num2str(length(idx_no_drizzle)), ' profiles'], 'Interpreter', 'latex', 'FontSize', ttl_fnt)

% N_c
subplot(1,3,3)
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
hold on
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

set(gcf, 'Position', [0 0 1200 600])

% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.213894399154591 0.144881944444444 0.189354393115943 0.771350694444444]);
hold(axes1,'on');

subplot(1,3,1)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.0153049178561151 0.812493412950304 0.171744791666667 0.114166666666667])



subplot(1,3,2)
hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)


subplot(1,3,3)
hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)

set(gcf, 'Position', [0 0 1200 600])

hold on; subplot(1,3,2)
title('ORACLES median profiles: drizzling vs non-drizzling', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)


%% -------------------------------------------------------------------------
%  Print summary statistics
%  -------------------------------------------------------------------------

fprintf('\n========== ORACLES Cloud Microphysics Summary ==========\n')
fprintf('Total profiles loaded:    %d\n', N_profiles)
fprintf('Non-drizzling profiles:   %d  (%.0f%%)\n', ...
    length(idx_no_drizzle), 100*length(idx_no_drizzle)/N_profiles)
fprintf('Drizzling profiles:       %d  (%.0f%%)\n', ...
    length(idx_drizzle),    100*length(idx_drizzle)/N_profiles)

fprintf('\n--- Adiabaticity (Wood 2005) ---\n')
fprintf('  Non-drizzling: median = %.2f  (IQR: %.2f – %.2f)\n', ...
    median(adiabaticity(~is_drizzle), 'omitnan'), ...
    quantile(adiabaticity(~is_drizzle), 0.25), ...
    quantile(adiabaticity(~is_drizzle), 0.75))
fprintf('  Drizzling:     median = %.2f  (IQR: %.2f – %.2f)\n', ...
    median(adiabaticity(is_drizzle), 'omitnan'), ...
    quantile(adiabaticity(is_drizzle), 0.25), ...
    quantile(adiabaticity(is_drizzle), 0.75))

fprintf('\n--- CTEI (kappa, lower-bound estimate) ---\n')
fprintf('  Non-drizzling: median kappa = %.2f  (kappa_c = %.2f)\n', ...
    median(kappa_CTEI(~is_drizzle), 'omitnan'), kappa_c)
fprintf('  Drizzling:     median kappa = %.2f\n', ...
    median(kappa_CTEI(is_drizzle), 'omitnan'))
fprintf('  NOTE: kappa estimated from moist-adiabat departure; ', ...
    'above-cloud data would improve this estimate.\n')

fprintf('\n--- Turbulent mixing regime (LWC vs tau_phase log-log slope) ---\n')
fprintf('  Non-drizzling: slope = %.2f\n', coeff_nd(1))
fprintf('  Drizzling:     slope = %.2f\n', coeff_d(1))
fprintf('  (homogeneous ≈ -1,  inhomogeneous ≈ 0)\n\n')


%% =========================================================================
%  Local functions
%  =========================================================================

function [re, lwc, Nc, alt, T_K, P_Pa, rv] = orientProfile(prof)
% ORIENTPROFILE  Orient a vertical profile to cloud-base-first order.
%
%   Outputs are all oriented so that index 1 corresponds to cloud base
%   and index end corresponds to cloud top, regardless of whether the
%   aircraft was ascending or descending.

dz_dt = mean(diff(prof.altitude) ./ diff(prof.time));

if dz_dt > 0
    % Ascending: data starts at cloud base — no flip needed
    re   = prof.re;
    lwc  = prof.lwc;
    Nc   = prof.total_Nc;    % cm⁻³
    alt  = prof.altitude;    % m MSL
    T_K  = prof.temp + 273.15;
    P_Pa = prof.pres * 100;  % mb → Pa
    rv   = prof.rv;          % µm
else
    % Descending: data starts at cloud top — flip to cloud-base-first
    re   = fliplr(prof.re);
    lwc  = fliplr(prof.lwc);
    Nc   = fliplr(prof.total_Nc);
    alt  = fliplr(prof.altitude);
    T_K  = fliplr(prof.temp + 273.15);
    P_Pa = fliplr(prof.pres * 100);
    rv   = fliplr(prof.rv);
end

% Convert row → column vectors for consistency
re   = re(:)';
lwc  = lwc(:)';
Nc   = Nc(:)';
alt  = alt(:)';
T_K  = T_K(:)';
P_Pa = P_Pa(:)';
rv   = rv(:)';

end
