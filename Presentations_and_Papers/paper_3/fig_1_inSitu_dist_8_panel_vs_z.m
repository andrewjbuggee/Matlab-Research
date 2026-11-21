%% Figure 1 - Combine VOCALS-REx and ORACLES - 10-90 as shaded bands
%
%
% Here we plot percentiles directly:
%   filled band  : 10th to 90th percentile
%   solid line   : median (50th percentile)
%



% -------------------------------------------------------------------------
% ------ VOCALS-REx: load, separate drizzle, compute percentiles ----------
% -------------------------------------------------------------------------
clear variables

tauC_min = 3;
tauC_max = 40;

% --- Effective radius cleaning ---------------------------------------------
% re_filter_mode chooses what happens to a level with r_e <= r_eff_min_um
% or r_e >= r_eff_max_um:
%
%   'levels'   : drop only that LEVEL and keep the rest of the profile.
%                This matches clean_insitu_profile() in
%                Hyperspectral_Cloud_Retrievals/Neural_Network/build_synthetic_profile_inputs_from_insitu_and_ERA5.m
%                (inputs.remove_levels_re_le_1um, inputs.remove_levels_re_ge_50um,
%                inputs.min_num_levels), whose output is what
%                06_synthetic_profile_generator.py trains on.
%   'profiles' : drop the WHOLE PROFILE if any of its levels is out of
%                range. Note that this also removes profiles with an r_e = 0
%                level (an empty bin, usually at a cloud edge).
%
% In both modes, levels with non-finite r_e, LWC or altitude are dropped,
% and a profile is skipped if fewer than min_num_levels levels remain.
re_filter_mode = 'levels';   % 'levels' or 'profiles'
r_eff_min_um   = 1;          % out of range if r_e <= this (microns)
r_eff_max_um   = 50;         % out of range if r_e >= this (microns)
min_num_levels = 3;          % skip a profile if fewer levels than this remain

if ~any(strcmp(re_filter_mode, {'levels', 'profiles'}))
    error([newline, 're_filter_mode must be ''levels'' or ''profiles''.', newline])
end

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

end

% Version from May 2026
% load([foldername_data,...
%     'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])

% Version from November 2026 with optical depth correction
load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.05_Nc-threshold_10_08-Nov-2026.mat'])


N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

n_bins     = 30;
bin_edges  = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

vertSeg_nd = cell(n_bins, 4);
vertSeg_d  = cell(n_bins, 4);

num_removed        = 0;   % profiles outside the tau_c window
num_re_profiles    = 0;   % profiles dropped for r_e ('profiles' mode only)
num_too_few_levels = 0;   % profiles with < min_num_levels after level cleaning
num_kept           = 0;
num_levels_removed = 0;   % total levels dropped by the r_e / finite-value cleaning

for nn = 1:N_profiles

    % --- optical depth filter (whole profile) ---
    if max(ensemble_profiles{nn}.tau) > tauC_max || max(ensemble_profiles{nn}.tau) < tauC_min
        num_removed = num_removed + 1;
        continue
    end

    % --- pick the effective radius the same way the build script does ---
    % Profiles whose 2DC data didn't conform lack the combined CDP + 2DC r_e;
    % fall back to the CDP-only r_e.
    if isfield(ensemble_profiles{nn}, 're')==true
        re_raw = ensemble_profiles{nn}.re;       % microns
    else
        re_raw = ensemble_profiles{nn}.re_CDP;   % microns
    end

    % --- r_e cleaning (see re_filter_mode at the top) ---
    % re_bad and keep_lvl are in the original sampling order of the profile
    re_bad = re_raw <= r_eff_min_um | re_raw >= r_eff_max_um;

    if strcmp(re_filter_mode, 'profiles') && any(re_bad)
        num_re_profiles = num_re_profiles + 1;
        continue
    end

    keep_lvl = ~re_bad & isfinite(re_raw) & ...
        isfinite(ensemble_profiles{nn}.lwc) & isfinite(ensemble_profiles{nn}.altitude);
    num_levels_removed = num_levels_removed + sum(~keep_lvl);

    if sum(keep_lvl) < min_num_levels
        num_too_few_levels = num_too_few_levels + 1;
        continue
    end
    num_kept = num_kept + 1;

    % Fit the gamma size distribution on all levels (the fit is per level,
    % so the result is subset below with the same mask)
    significance_lvl = 0.1;
    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nc, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);
    v_eff_all = 1 ./ (3 + gammaFit.alpha);   % column vector, one per level

    % keep only the cleaned levels
    z_m       = ensemble_profiles{nn}.altitude(keep_lvl);   % meters
    time_s    = ensemble_profiles{nn}.time(keep_lvl);       % seconds
    re_bin    = re_raw(keep_lvl);
    lwc_bin   = ensemble_profiles{nn}.lwc(keep_lvl);
    Nc_bin    = ensemble_profiles{nn}.total_Nc(keep_lvl);
    v_eff_bin = v_eff_all(keep_lvl);

    % Normalize altitude over the levels that remain, so 0 and 1 are the
    % lowest and highest kept levels. This matches the generator, which
    % maps each cleaned profile onto its own normalized-altitude grid.
    z_norm = (z_m - min(z_m)) ./ (max(z_m) - min(z_m));

    % orient so the arrays run from cloud top to cloud base
    dz_dt = mean(diff(z_m) ./ diff(time_s));
    if dz_dt > 0
        re_bin    = fliplr(re_bin);
        lwc_bin   = fliplr(lwc_bin);
        Nc_bin    = fliplr(Nc_bin);
        v_eff_bin = flipud(v_eff_bin);
        z_norm    = fliplr(z_norm);
    end

    % r_e <= 1 um is already removed, so only LWC and Nc can still be zero
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0)   = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = z_norm >= bin_edges(bb) & z_norm <= bin_edges(bb+1);
        else
            idx_seg = z_norm >  bin_edges(bb) & z_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end

% Percentiles per bin: [p10, p50, p90]
%   column 1 (p10) and column 3 (p90) : edges of the shaded band
%   column 2 (p50)                    : median, plotted as the solid line
pct = [10 50 90];
re_p_nd   = zeros(n_bins, length(pct));  re_p_d   = zeros(n_bins, length(pct));
lwc_p_nd  = zeros(n_bins, length(pct));  lwc_p_d  = zeros(n_bins, length(pct));
Nc_p_nd   = zeros(n_bins, length(pct));  Nc_p_d   = zeros(n_bins, length(pct));
vEff_p_nd = zeros(n_bins, length(pct));  vEff_p_d = zeros(n_bins, length(pct));

for bb = 1:n_bins
    re_p_nd(bb,:)   = prctile(vertSeg_nd{bb,1}, pct);
    lwc_p_nd(bb,:)  = prctile(vertSeg_nd{bb,2}, pct);
    Nc_p_nd(bb,:)   = prctile(vertSeg_nd{bb,3}, pct);
    vEff_p_nd(bb,:) = prctile(vertSeg_nd{bb,4}, pct);

    re_p_d(bb,:)    = prctile(vertSeg_d{bb,1}, pct);
    lwc_p_d(bb,:)   = prctile(vertSeg_d{bb,2}, pct);
    Nc_p_d(bb,:)    = prctile(vertSeg_d{bb,3}, pct);
    vEff_p_d(bb,:)  = prctile(vertSeg_d{bb,4}, pct);
end

disp(['r_e filter mode            : ', re_filter_mode])
disp(['# profiles kept            : ', num2str(num_kept)])
disp(['# removed, tau_c window    : ', num2str(num_removed)])
disp(['# removed, r_e (profiles)  : ', num2str(num_re_profiles)])
disp(['# removed, too few levels  : ', num2str(num_too_few_levels)])
disp(['# levels removed, r_e/NaN  : ', num2str(num_levels_removed)])

% -------------------------------------------------------------------------
% ------ Plot VOCALS-REx panels (a)-(d) -----------------------------------
% -------------------------------------------------------------------------
figure1 = figure;

plt_clr_1   = mySavedColors(68,'fixed');
plt_clr_2   = mySavedColors(67,'fixed');

fnt_sz      = 18;
ttl_fnt     = 20;
alpha_inner = 0.30;
ln_w        = 1.5;
pnl_fnt     = 21;

re_lim_min = 3;   % microns
re_lim_max = 30;      % microns

y = [bin_center; flipud(bin_center)];

% --- (a) r_e ---
subplot(2,4,1)
hold on
% plot 10-90% data range as shaded region
fill([re_p_nd(:,1); flipud(re_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([re_p_d(:,1);  flipud(re_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([re_lim_min, re_lim_max])

annotation(figure1,'textbox',[0.130138 0.836 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (b) v_eff ---
subplot(2,4,2)
hold on
% plot 10-90% data range as shaded region
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

set(gca, 'XScale', 'log')
grid on; grid minor
xlim([0.005, 0.175])
% MATLAB's automatic log ticks only place 10^-2 on this ~1.5-decade range,
% even though 0.1 is inside the limits, so set both decades explicitly
xticks([0.01, 0.1])

annotation(figure1,'textbox',[0.3351315023418542 0.836 0.0373 0.0666],'String',{'(b)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (c) LWC ---
subplot(2,4,3)
hold on
% plot 10-90% data range as shaded region
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
% title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_min)],...
%     'Interpreter', 'latex','FontSize', ttl_fnt)
title(['VOCALS-REx median profiles - N = ', num2str(num_kept)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
annotation(figure1,'textbox',[0.540879493143727 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');
xlim([0, 0.725])

% --- (d) N_c ---
subplot(2,4,4)
hold on
% plot 10-90% data range as shaded region
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlim([0, 325])
annotation(figure1,'textbox',[0.746454846227638 0.836 0.0373 0.0666],'String','(d)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

legend('', 'Non-drizzling','', 'Drizzling', ...
    'Location','best','Interpreter','latex', ...
    'FontSize', fnt_sz-2, 'Color','white','TextColor','black', ...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])


% -------------------------------------------------------------------------
% ------ ORACLES: load, separate drizzle, compute percentiles -------------
% -------------------------------------------------------------------------
clearvars -except figure1 bin_center n_bins bin_edges fnt_sz ttl_fnt ...
    plt_clr_1 plt_clr_2 alpha_inner ln_w pnl_fnt re_lim_max re_lim_min...
    tauC_min tauC_max re_filter_mode r_eff_min_um r_eff_max_um min_num_levels

tauC_min = 3;

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

end

% % Version from May 2026
% load([foldername_data,...
%     'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

% Version from November 2026 with optical depth correction
load([foldername_data,...
    'ensemble_profiles_from_33_files_LWC-threshold_0.05_Nc-threshold_10_08-Nov-2026.mat'])


N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

% --- Outliers excluded from the synthetic profile generator ---------------
% These are the same two ORACLES profiles listed in OUTLIERS in
% lasp-CU-paper-3/06_synthetic_profile_generator.py (lines 103-151). There,
% profile_idx is the 1-based position in ensemble_profiles of the
% 08-Nov-2026 ORACLES file loaded above (it is 'kk' in
% build_synthetic_profile_inputs_from_insitu_and_ERA5.m), so the indices
% carry over directly.
%
% Both are from 30 Aug 2017: median LWC of 1.9 and 2.3 g/m^3 in layers only
% 122 and 125 m deep, i.e. f_ad = LWP_obs/LWP_ad of 9.3 and 10.4 (Wood 2005
% definition). Every other profile in the generator's training set has
% f_ad <= 2.75. See the generator header for the full justification.
% Empty these two arrays to restore them.
outlier_idx  = [111, 112];
outlier_date = {'2017-08-30', '2017-08-30'};

% Guard: make sure the indices still point at the intended profiles. If
% the ensemble file is rebuilt the indices can shift silently.
for oo = 1:length(outlier_idx)
    date_oo = char(ensemble_profiles{outlier_idx(oo)}.dateOfFlight, 'yyyy-MM-dd');
    if ~strcmp(date_oo, outlier_date{oo})
        error([newline, 'ORACLES outlier index ', num2str(outlier_idx(oo)), ...
            ' has flight date ', date_oo, ', expected ', outlier_date{oo}, '.', ...
            newline, 'The ensemble file may have changed; update outlier_idx.', newline])
    end
end

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

vertSeg_nd = cell(n_bins, 4);
vertSeg_d  = cell(n_bins, 4);

num_removed        = 0;   % profiles outside the tau_c window
num_re_profiles    = 0;   % profiles dropped for r_e ('profiles' mode only)
num_too_few_levels = 0;   % profiles with < min_num_levels after level cleaning
num_kept           = 0;
num_levels_removed = 0;   % total levels dropped by the r_e / finite-value cleaning
num_outliers       = 0;

for nn = 1:N_profiles

    % skip the outliers first so they aren't also counted in num_removed
    if ismember(nn, outlier_idx)
        num_outliers = num_outliers + 1;
        continue
    end

    % --- optical depth filter (whole profile) ---
    if max(ensemble_profiles{nn}.tau) > tauC_max || max(ensemble_profiles{nn}.tau) < tauC_min
        num_removed = num_removed + 1;
        continue
    end

    re_raw = ensemble_profiles{nn}.re;   % microns (CAS + 2DS/HVPS)

    % --- r_e cleaning (see re_filter_mode in the VOCALS-REx section) ---
    % re_bad and keep_lvl are in the original sampling order of the profile
    re_bad = re_raw <= r_eff_min_um | re_raw >= r_eff_max_um;

    if strcmp(re_filter_mode, 'profiles') && any(re_bad)
        num_re_profiles = num_re_profiles + 1;
        continue
    end

    keep_lvl = ~re_bad & isfinite(re_raw) & ...
        isfinite(ensemble_profiles{nn}.lwc) & isfinite(ensemble_profiles{nn}.altitude);
    num_levels_removed = num_levels_removed + sum(~keep_lvl);

    if sum(keep_lvl) < min_num_levels
        num_too_few_levels = num_too_few_levels + 1;
        continue
    end
    num_kept = num_kept + 1;

    % Fit the gamma size distribution on all levels (the fit is per level,
    % so the result is subset below with the same mask)
    significance_lvl = 0.1;
    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nd, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);
    v_eff_all = 1 ./ (3 + gammaFit.alpha);   % column vector, one per level

    % keep only the cleaned levels
    z_m       = ensemble_profiles{nn}.altitude(keep_lvl);   % meters
    time_s    = ensemble_profiles{nn}.time(keep_lvl);       % seconds
    re_bin    = re_raw(keep_lvl);
    lwc_bin   = ensemble_profiles{nn}.lwc(keep_lvl);
    Nc_bin    = ensemble_profiles{nn}.total_Nc(keep_lvl);
    v_eff_bin = v_eff_all(keep_lvl);

    % normalize altitude over the levels that remain (see VOCALS-REx section)
    z_norm = (z_m - min(z_m)) ./ (max(z_m) - min(z_m));

    % orient so the arrays run from cloud top to cloud base
    dz_dt = mean(diff(z_m) ./ diff(time_s));
    if dz_dt > 0
        re_bin    = fliplr(re_bin);
        lwc_bin   = fliplr(lwc_bin);
        Nc_bin    = fliplr(Nc_bin);
        v_eff_bin = flipud(v_eff_bin);
        z_norm    = fliplr(z_norm);
    end

    % r_e <= 1 um is already removed, so only LWC and Nc can still be zero
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0)   = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = z_norm >= bin_edges(bb) & z_norm <= bin_edges(bb+1);
        else
            idx_seg = z_norm >  bin_edges(bb) & z_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end

% Percentiles per bin: [p10, p50, p90] (same layout as the VOCALS-REx block)
pct = [10 50 90];
re_p_nd   = zeros(n_bins, length(pct));  re_p_d   = zeros(n_bins, length(pct));
lwc_p_nd  = zeros(n_bins, length(pct));  lwc_p_d  = zeros(n_bins, length(pct));
Nc_p_nd   = zeros(n_bins, length(pct));  Nc_p_d   = zeros(n_bins, length(pct));
vEff_p_nd = zeros(n_bins, length(pct));  vEff_p_d = zeros(n_bins, length(pct));

for bb = 1:n_bins
    re_p_nd(bb,:)   = prctile(vertSeg_nd{bb,1}, pct);
    lwc_p_nd(bb,:)  = prctile(vertSeg_nd{bb,2}, pct);
    Nc_p_nd(bb,:)   = prctile(vertSeg_nd{bb,3}, pct);
    vEff_p_nd(bb,:) = prctile(vertSeg_nd{bb,4}, pct);

    re_p_d(bb,:)    = prctile(vertSeg_d{bb,1}, pct);
    lwc_p_d(bb,:)   = prctile(vertSeg_d{bb,2}, pct);
    Nc_p_d(bb,:)    = prctile(vertSeg_d{bb,3}, pct);
    vEff_p_d(bb,:)  = prctile(vertSeg_d{bb,4}, pct);
end

disp(['r_e filter mode            : ', re_filter_mode])
disp(['# profiles kept            : ', num2str(num_kept)])
disp(['# removed, tau_c window    : ', num2str(num_removed)])
disp(['# removed, r_e (profiles)  : ', num2str(num_re_profiles)])
disp(['# removed, too few levels  : ', num2str(num_too_few_levels)])
disp(['# levels removed, r_e/NaN  : ', num2str(num_levels_removed)])
disp(['# outliers removed         : ', num2str(num_outliers)])

y = [bin_center; flipud(bin_center)];

% -------------------------------------------------------------------------
% ------ Plot ORACLES panels (e)-(h) --------------------------------------
% -------------------------------------------------------------------------

% --- (e) r_e ---
subplot(2,4,5)
hold on
% plot 10-90% data range as shaded region
fill([re_p_nd(:,1); flipud(re_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([re_p_d(:,1);  flipud(re_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([re_lim_min, re_lim_max])
annotation(figure1,'textbox',[0.132465746741153 0.367111111111112 0.0373 0.0666],'String',{'(e)'},...
    'Interpreter','latex','FontSize', pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (f) v_eff ---
subplot(2,4,6)
hold on
% plot 10-90% data range as shaded region
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

set(gca, 'XScale', 'log')
grid on; grid minor
xlabel('$\langle \nu_{e}(z) \rangle$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.005, 0.175])
% set both decades explicitly (see panel (b))
xticks([0.01, 0.1])
annotation(figure1,'textbox',[0.343867473506009 0.367111111111112 0.0373 0.0666],'String',{'(f)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (g) LWC ---
subplot(2,4,7)
hold on
% plot 10-90% data range as shaded region
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_min)],...
%     'Interpreter', 'latex','FontSize', ttl_fnt)
title(['ORACLES median profiles - N = ', num2str(num_kept)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.725])
annotation(figure1,'textbox',[0.542078635404322 0.361555555555556 0.0372999999999999 0.0666],'String','(g)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (h) N_c ---
subplot(2,4,8)
hold on
% plot 10-90% data range as shaded region
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 325])
annotation(figure1,'textbox',[0.745932819705434 0.367111111111112 0.0373 0.0666],'String','(h)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% IEEE-column sizing
w = 7.16; h = 3;
figure1.Units    = 'inches';
figure1.Position = [1, 1, 2.75*w, 2.75*h];


% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle and non drizzle',...
    ' - 8 panels with 10 and 90 percentiles - function of z.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle',...
    'and non drizzle - 8 panels with 10 and 90 percentiles - function of z.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------
