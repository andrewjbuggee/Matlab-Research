%% Wood (2005) adiabaticity: look at how cloud base and Gamma_ad are chosen
%
% WHAT THIS SCRIPT IS FOR
%   compute_adiabaticity_Wood2005 estimates f_ad = LWP_obs / LWP_ad for every
%   in-situ profile, with LWP_ad = (1/2) Gamma_ad h^2 and h = z_top - z_cb. Two
%   choices decide the answer, and this script draws both so they can be checked
%   by eye:
%
%   CLOUD BASE z_cb. A straight line is fit to LWC(z) over the bottom third of
%   each profile, and its z-intercept (LWC = 0) is the cloud base, but never above
%   the lowest kept sample: if the line reaches zero above it, the lowest sample is
%   the base. If the fitted slope is not physical (<= 0.25 Gamma_ad or > 3 Gamma_ad,
%   or fewer than 3 samples) the fallback puts the base one adiabatic distance below
%   the lowest sample:
%       z_cb = z_lowest - LWC_threshold / Gamma_ad
%
%   GAMMA_AD. One value per profile, at the mean in-cloud temperature and pressure
%   (Wood 2005, p. 3019), from the moist-adiabatic relations written out in
%   Figure 6.
%
% FIGURES
%   1. Why the slope has limits: the fitted slopes, and how far below the cloud
%      the intercept would land without a lower limit.
%   2. Every profile's bottom-third fit (pages of panels).
%   3. The fallback cases, each labeled with why its fit was rejected.
%   4. Profiles whose line reaches LWC = 0 ABOVE the lowest sample, which therefore
%      use the lowest sample as the base: z_cb and f_ad against the intercept.
%   5. The LWC profiles behind Figure 4.
%   6. Mean in-cloud T and p, the Gamma_ad computed from them, and the equations.
%   7. Climbs against descents: in-cloud T, p, Gamma_ad and the pieces of f_ad.
%   8. Where in the cloud climbs and descents differ: LWC / q_ad against z*.
%   9. Two tests of where a climb-descent difference could come from: neighboring
%      climb-descent pairs, and a time shift between the LWC and altitude records.
%
%   Nothing is refit here. The profiles come from the same profile finders and
%   thresholds as the ensemble scripts, and every line, base and f_ad is read
%   from the .adiabaticity structures those finders store, so the figures show
%   exactly what the pipeline uses.
%
% SEE ALSO:
%   compute_adiabaticity_Wood2005, find_verticalProfiles_VOCALS_REx_ver2,
%   find_verticalProfiles_ORACLES, plot_adiabaticity_base_comparison
%
% By Andrew John Buggee
%%

clear variables


%% Settings

% Profile thresholds: the same values as ensemble_vertical_statistics.m and
% ensemble_vertical_statistics_ORACLES.m, so the same profiles are found
inputs.LWC_threshold = 0.05;        % g/m^3
inputs.Nc_threshold = 10;           % cm^-3
inputs.stop_at_max_LWC = false;

% ORACLES is needed for the campaign comparison in Figures 4-6. Turning it off
% skips about 25 s of reading.
inputs.include_oracles = true;

% Which campaign gets the page-by-page fits of Figure 2 and the fallback panels
% of Figure 3: 'VOCALS-REx' (92 profiles, 4 pages) or 'ORACLES' (239, 10 pages)
inputs.fit_pages_campaign = 'VOCALS-REx';
inputs.panels_per_page = [4, 6];    % [rows, columns]

% '' shows every figure on screen. A folder path saves PNGs there instead; the
% figures are then created hidden and closed after saving, because R2024b can hang
% while exporting a visible figure.
inputs.save_figures_to = '';

% Climb-descent comparison (Figures 7-9)
inputs.rng_seed = 42;                       % seed of the bootstrap confidence intervals
inputs.n_bootstrap = 5000;                  % bootstrap resamples
inputs.pair_max_separation_s = 900;         % climb and descent closer than this count as a pair
inputs.lwc_time_shifts_s = -6:6;            % s; shifts of the LWC record against altitude (Figure 9)


% ------------------------------------------------------------------------
% Data folders
% ------------------------------------------------------------------------
which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    folder_vocals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
    folder_oracles = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/combined_microphysics/'];

elseif strcmp(which_computer, 'andrewbuggee')

    folder_vocals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];
    folder_oracles = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/ORACLES_data/combined_microphysics/'];

end


% ------------------------------------------------------------------------
% Colors. One role per color across every figure in this script. The trio was
% checked with the dataviz palette validator (all pairs pass for color-vision
% deficiency); aqua and yellow are below 3:1 contrast on white, so they always
% come with a legend entry or a printed table.
% ------------------------------------------------------------------------
style.fit        = [0.922, 0.408, 0.204];   % orange #eb6834  line fit and its intercept base
style.fallback   = [0.106, 0.686, 0.478];   % aqua   #1baf7a  fallback base
style.lowest     = [0.165, 0.471, 0.839];   % blue   #2a78d6  lowest sample taken as the base
style.campaign   = [0.384, 0.314, 0.839;    % violet #6250d6  VOCALS-REx
                    0.929, 0.631, 0.000];   % yellow #eda100  ORACLES
style.gray_light = [0.72, 0.72, 0.72];      % samples above the fitted layer
style.gray_mid   = [0.55, 0.55, 0.55];      % reference rules
style.gray_dark  = [0.20, 0.20, 0.20];      % adiabatic profile


%% Find the profiles and read their adiabaticity (same finders as the ensemble scripts)

% The VOCALS-REx ensemble script also sorts profiles for drizzle, but it is set to
% keep the drizzling ones, so every profile the finder returns is in its ensemble.
records = {};

files_vocals = dir([folder_vocals, 'RF*.nc']);

for ff = 1:numel(files_vocals)

    vocalsRex = readVocalsRex([folder_vocals, files_vocals(ff).name]);
    vert_profs = find_verticalProfiles_VOCALS_REx_ver2(vocalsRex, inputs.LWC_threshold, ...
        inputs.stop_at_max_LWC, inputs.Nc_threshold, which_computer);

    for nn = 1:numel(vert_profs)
        % Label: research flight and the profile's position in the 92-profile ensemble
        label = sprintf('%s #%d', files_vocals(ff).name(1:4), numel(records) + 1);
        rec = profile_record('VOCALS-REx', label, files_vocals(ff).name(1:4), vert_profs(nn));
        % f_ad again with the LWC record shifted in time against altitude (Figure 9)
        rec.f_ad_by_shift = f_ad_with_time_shift(vocalsRex.time, vocalsRex.altitude, vocalsRex.lwc, ...
            vocalsRex.temp, vocalsRex.pres, vert_profs(nn).time, inputs.LWC_threshold, inputs.lwc_time_shifts_s);
        records{end+1} = rec; %#ok<SAGROW>
    end

end

n_vocals = numel(records);
disp([newline, 'VOCALS-REx: ', num2str(n_vocals), ' profiles'])

if inputs.include_oracles

    files_oracles = dir([folder_oracles, 'Microphysics_P*.nc']);

    for ff = 1:numel(files_oracles)

        oracles = readORACLES([folder_oracles, files_oracles(ff).name]);
        vert_profs = find_verticalProfiles_ORACLES(oracles, inputs.LWC_threshold, ...
            inputs.stop_at_max_LWC, inputs.Nc_threshold, which_computer);

        % Some flights have no profile that passes the thresholds
        if isempty(vert_profs)
            continue
        end

        flight_id = char(datetime(vert_profs(1).dateOfFlight, 'Format', 'dd-MMM-yy'));

        for nn = 1:numel(vert_profs)
            label = sprintf('%s #%d', flight_id, numel(records) - n_vocals + 1);
            rec = profile_record('ORACLES', label, flight_id, vert_profs(nn));
            rec.f_ad_by_shift = f_ad_with_time_shift(oracles.time, oracles.altitude, oracles.lwc, ...
                oracles.temp, oracles.pres, vert_profs(nn).time, inputs.LWC_threshold, inputs.lwc_time_shifts_s);
            records{end+1} = rec; %#ok<SAGROW>
        end

    end

    disp(['ORACLES: ', num2str(numel(records) - n_vocals), ' profiles'])

end

profiles = [records{:}];


%% One row per profile: everything the figures below use

n_prof = numel(profiles);

campaign            = string({profiles.campaign})';
label               = string({profiles.label})';
is_ascending        = [profiles.is_ascending]';
T_mean_C            = arrayfun(@(p) p.a.T_mean_C, profiles)';
P_mean_hPa          = arrayfun(@(p) p.a.P_mean_hPa, profiles)';
Gamma_ad_g_m3_per_km = arrayfun(@(p) p.a.Gamma_ad_g_m3_per_km, profiles)';
slope_over_gamma_ad = arrayfun(@(p) p.a.fit_slope_over_gamma_ad, profiles)';
n_fit_samples       = arrayfun(@(p) p.a.n_fit_samples, profiles)';
base_source         = string(arrayfun(@(p) p.a.base_source, profiles, 'UniformOutput', false))';
z_cb_above_lowest_m = arrayfun(@(p) -p.a.dz_base_m, profiles)';
h_m                 = arrayfun(@(p) p.a.h_m, profiles)';
wedge_fraction      = arrayfun(@(p) p.a.lwp_wedge_g_m2 ./ p.a.lwp_obs_g_m2, profiles)';
f_ad                = arrayfun(@(p) p.a.f_ad_lwp, profiles)';
f_ad_lowest_sample  = arrayfun(@(p) p.a_lowest.f_ad_lwp, profiles)';
lwp_obs_g_m2        = arrayfun(@(p) p.a.lwp_obs_g_m2, profiles)';
abs_w_m_s           = abs([profiles.w_m_s])';
flight              = string({profiles.flight})';
t_start_s           = [profiles.t_start_s]';
z_top_m             = arrayfun(@(p) p.a.z_top_m, profiles)';
z_lowest_m          = arrayfun(@(p) p.a.z_lowest_sample_m, profiles)';

% Where the fitted line itself reaches LWC = 0, whether or not it was used [m above the
% lowest sample]. Only defined for positive slopes.
line_zero_above_lowest_m = arrayfun(@(p) -p.a.lwc_fit_at_lowest_g_m3 ./ (p.a.fit_slope_g_m3_per_km ./ 1e3), profiles)';
line_zero_above_lowest_m(~(slope_over_gamma_ad > 0)) = NaN;

fallback_reason = strings(n_prof, 1);
for nn = 1:n_prof
    fallback_reason(nn) = reason_for_fallback(profiles(nn).a);
end

summary = table(campaign, label, is_ascending, T_mean_C, P_mean_hPa, Gamma_ad_g_m3_per_km, ...
    slope_over_gamma_ad, n_fit_samples, base_source, fallback_reason, z_cb_above_lowest_m, ...
    line_zero_above_lowest_m, h_m, wedge_fraction, f_ad, f_ad_lowest_sample, lwp_obs_g_m2, ...
    abs_w_m_s, flight, t_start_s, z_top_m, z_lowest_m);

is_fallback = summary.base_source == "fallback";
is_lowest_base = summary.base_source == "lowest sample, intercept above it";
is_vocals   = summary.campaign == "VOCALS-REx";
is_oracles  = summary.campaign == "ORACLES";
has_gamma   = isfinite(summary.Gamma_ad_g_m3_per_km);
campaign_names = ["VOCALS-REx", "ORACLES"];
campaign_short = ["VOCALS", "ORACLES"];

for cc = ["VOCALS-REx", "ORACLES"]
    k = summary.campaign == cc;
    if any(k)
        fprintf(['\n%s: %d profiles, %d fallbacks, %d whose line reaches zero above the lowest sample ', ...
            '(that sample is the base), %d without temperature.\n', ...
            '  Unsampled wedge as a fraction of LWP_obs: median %.3f, max %.3f\n'], ...
            cc, sum(k), sum(k & is_fallback), sum(k & is_lowest_base), ...
            sum(k & ~has_gamma), median(summary.wedge_fraction(k), 'omitnan'), max(summary.wedge_fraction(k)));
    end
end


%% Figure 1: why the fitted slope has a lower limit above zero

fig = make_figure(inputs, 'Wood adiabaticity: slope limits', [40, 40, 1300, 520]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
limits = profiles(1).a.slope_limits_over_gamma_ad;

% ---- (a) Distribution of the fitted slope ---------------------------------
ax = nexttile(tl); hold(ax, 'on')

edges = -1.5:0.1:4;
y_max = 0;
h = gobjects(0); names = {};
for cc = 1:2
    k = has_gamma & summary.campaign == campaign_names(cc);
    if ~any(k), continue, end
    counts = histcounts(min(max(summary.slope_over_gamma_ad(k), edges(1)), edges(end)), edges);
    h(end+1) = stairs(ax, edges, [counts, counts(end)], 'Color', style.campaign(cc, :), 'LineWidth', 2); %#ok<SAGROW>
    names{end+1} = sprintf('%s: %d of %d fall back', campaign_names(cc), ...
        sum(k & is_fallback), sum(k)); %#ok<SAGROW>
    y_max = max(y_max, max(counts));
end

% Shade the slopes that send a profile to the fallback
patch(ax, [edges(1), limits(1), limits(1), edges(1)], [0, 0, 1.1, 1.1] .* y_max, style.gray_mid, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
patch(ax, [limits(2), edges(end), edges(end), limits(2)], [0, 0, 1.1, 1.1] .* y_max, style.gray_mid, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
xline(ax, [0, 1], 'Color', style.gray_mid, 'LineWidth', 0.5, 'HandleVisibility', 'off')
text(ax, limits(1) - 0.05, 1.05 * y_max, sprintf('fallback\n(slope \\leq %.2f \\Gamma_{ad})', limits(1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9, 'Color', style.gray_dark)
text(ax, limits(2) + 0.05, 1.05 * y_max, sprintf('fallback\n(slope > %g \\Gamma_{ad})', limits(2)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 9, 'Color', style.gray_dark)
text(ax, 1.03, 0.62 * y_max, 'adiabatic', 'FontSize', 9, 'Color', style.gray_dark)

legend(ax, h, names, 'Location', 'east', 'Box', 'off')
xlim(ax, edges([1, end])); ylim(ax, [0, 1.1 * y_max])
xlabel(ax, 'slope of the line fit to the bottom third  /  \Gamma_{ad}')
ylabel(ax, 'profiles')
title(ax, '(a) Fitted slope, in units of each profile''s \Gamma_{ad} (end bins hold everything beyond)', ...
    'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1;

% ---- (b) How far below the cloud the line reaches zero ---------------------
ax = nexttile(tl); hold(ax, 'on')

h = gobjects(0); names = {};
for cc = 1:2
    k = has_gamma & summary.campaign == campaign_names(cc) & summary.line_zero_above_lowest_m < 0;
    if ~any(k), continue, end
    h(end+1) = scatter(ax, summary.slope_over_gamma_ad(k), -summary.line_zero_above_lowest_m(k), 22, ...
        style.campaign(cc, :), 'filled', 'MarkerFaceAlpha', 0.75); %#ok<SAGROW>
    names{end+1} = char(campaign_names(cc)); %#ok<SAGROW>
end

k = has_gamma & is_fallback & summary.line_zero_above_lowest_m < 0;
h(end+1) = scatter(ax, summary.slope_over_gamma_ad(k), -summary.line_zero_above_lowest_m(k), 70, ...
    style.fallback, 'x', 'LineWidth', 1.5);
names{end+1} = 'fallback used instead';

xline(ax, limits, 'Color', style.gray_mid, 'LineWidth', 1, 'HandleVisibility', 'off')
set(ax, 'XScale', 'log', 'YScale', 'log')
legend(ax, h, names, 'Location', 'northeast', 'Box', 'off')
xlabel(ax, 'slope of the line fit to the bottom third  /  \Gamma_{ad}')
ylabel(ax, 'depth below the lowest sample where the line reaches LWC = 0 (m)')

kept = has_gamma & ~is_fallback & summary.z_cb_above_lowest_m < 0;
title(ax, sprintf(['(b) Depth = LWC_{line}(z_{lowest}) / slope, which grows without bound as the ', ...
    'slope \\rightarrow 0.\nDeepest base without a lower limit: %.0f m; with it: %.0f m'], ...
    max(-summary.line_zero_above_lowest_m(has_gamma)), max(-summary.z_cb_above_lowest_m(kept))), ...
    'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1;

finish_figure(fig, inputs, 'wood_adiab_fig1_slope_limits')


%% Figure 2: every profile's bottom-third fit, page by page

k_pages = find(summary.campaign == inputs.fit_pages_campaign);
per_page = prod(inputs.panels_per_page);
n_pages = ceil(numel(k_pages) / per_page);

for pp = 1:n_pages

    idx = k_pages((pp - 1) * per_page + 1 : min(pp * per_page, numel(k_pages)));

    fig = make_figure(inputs, sprintf('Bottom-third fits, page %d', pp), [30, 30, 1700, 1050]);
    tl = tiledlayout(fig, inputs.panels_per_page(1), inputs.panels_per_page(2), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for jj = 1:numel(idx)
        ax = nexttile(tl);
        draw_profile_panel(ax, profiles(idx(jj)), style, 'fit');
        if jj == 1
            ax_legend = ax;
        end
    end

    add_panel_legend(ax_legend, style, any(is_fallback(idx)), any(is_lowest_base(idx)));
    xlabel(tl, 'LWC (g m^{-3})')
    ylabel(tl, 'height above the lowest kept sample (m)')
    title(tl, sprintf(['%s: line fit to LWC(z) over the bottom third; its intercept is z_{cb}. ', ...
        'Page %d of %d'], inputs.fit_pages_campaign, pp, n_pages))

    finish_figure(fig, inputs, sprintf('wood_adiab_fig2_fits_page%02d', pp))

end


%% Figure 3: the fallback cases

idx = find(summary.campaign == inputs.fit_pages_campaign & is_fallback);

disp([newline, 'Fallback cases, ', inputs.fit_pages_campaign, ':'])
disp(summary(idx, {'label', 'is_ascending', 'slope_over_gamma_ad', 'n_fit_samples', 'fallback_reason', ...
    'line_zero_above_lowest_m', 'z_cb_above_lowest_m', 'h_m', 'wedge_fraction', 'f_ad'}))

if ~isempty(idx)

    n_cols = min(5, numel(idx));
    n_rows = ceil(numel(idx) / n_cols);

    fig = make_figure(inputs, 'Fallback cases', [30, 30, 360 * n_cols, 330 * n_rows + 120]);
    tl = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for jj = 1:numel(idx)
        ax = nexttile(tl);
        draw_profile_panel(ax, profiles(idx(jj)), style, 'fallback');
        if jj == 1
            ax_legend = ax;
        end
    end

    add_panel_legend(ax_legend, style, true, false);
    xlabel(tl, 'LWC (g m^{-3})')
    ylabel(tl, 'height above the lowest kept sample (m)')
    title(tl, sprintf(['%s: %d profiles whose fitted slope was rejected. The dotted line is the rejected fit; ', ...
        'the base is put one adiabatic distance below the lowest sample instead'], ...
        inputs.fit_pages_campaign, numel(idx)))

    finish_figure(fig, inputs, 'wood_adiab_fig3_fallbacks')

end


%% Figure 4: line reaches LWC = 0 above the lowest sample, so the lowest sample is the base
% compute_adiabaticity_Wood2005 never puts the base above a measured sample. For these
% profiles the fitted line reaches LWC = 0 a few metres above the lowest kept sample,
% so the lowest sample is the base and every sample counts toward LWP_obs. To show
% how much that choice matters, the alternative is computed here from the same line
% and Gamma_ad: the intercept as the base, with the samples below it left out of
% LWP_obs (f_ad_with_intercept_base below).

idx = find(is_lowest_base);

h_intercept_m = nan(numel(idx), 1);
lwp_obs_intercept = nan(numel(idx), 1);
lwp_ad_intercept = nan(numel(idx), 1);
f_ad_intercept = nan(numel(idx), 1);
for jj = 1:numel(idx)
    [f_ad_intercept(jj), h_intercept_m(jj), lwp_obs_intercept(jj), lwp_ad_intercept(jj)] = ...
        f_ad_with_intercept_base(profiles(idx(jj)));
end

lwp_obs_used = arrayfun(@(p) p.a.lwp_obs_g_m2, profiles(idx))';
lwp_ad_used  = arrayfun(@(p) p.a.lwp_ad_g_m2, profiles(idx))';

lowest_base = table(summary.campaign(idx), summary.label(idx), summary.is_ascending(idx), ...
    summary.line_zero_above_lowest_m(idx), summary.h_m(idx), h_intercept_m, ...
    lwp_obs_used, lwp_obs_intercept, lwp_ad_used, lwp_ad_intercept, summary.f_ad(idx), f_ad_intercept, ...
    'VariableNames', {'campaign', 'label', 'is_ascending', 'line_zero_above_lowest_m', 'h_used_m', ...
    'h_intercept_m', 'lwp_obs_used', 'lwp_obs_intercept', 'lwp_ad_used', 'lwp_ad_intercept', ...
    'f_ad_used', 'f_ad_intercept'});

disp([newline, 'Lowest sample used as the base because the line reaches LWC = 0 above it (LWP in g/m^2):'])
disp(lowest_base)

if ~isempty(idx)

    fig = make_figure(inputs, 'Lowest sample as the base', [40, 40, 1300, 120 + 34 * numel(idx)]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    y_pos = numel(idx):-1:1;
    row_labels = flip(compose('%s %s (%s)', lowest_base.campaign, lowest_base.label, ...
        string(direction_word(lowest_base.is_ascending))));

    % ---- (a) cloud base relative to the lowest sample ----
    ax = nexttile(tl); hold(ax, 'on')
    for jj = 1:numel(idx)
        plot(ax, [0, lowest_base.line_zero_above_lowest_m(jj)], [y_pos(jj), y_pos(jj)], '-', ...
            'Color', style.gray_light, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
    h1 = plot(ax, zeros(size(y_pos)), y_pos, 's', 'MarkerSize', 8, ...
        'MarkerFaceColor', style.lowest, 'MarkerEdgeColor', style.lowest);
    h2 = plot(ax, lowest_base.line_zero_above_lowest_m, y_pos, '^', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', style.fit, 'LineWidth', 1.5);
    yticks(ax, 1:numel(idx)); yticklabels(ax, row_labels)
    ylim(ax, [0.4, numel(idx) + 0.6])
    xlabel(ax, 'height above the lowest kept sample (m)')
    legend(ax, [h1, h2], {'cloud base used: lowest sample', 'where the line reaches LWC = 0 (not used)'}, ...
        'Location', 'southeast', 'Box', 'off')
    title(ax, '(a) Cloud base', 'FontWeight', 'normal')
    grid(ax, 'on'); ax.GridAlpha = 0.1; ax.YGrid = 'off';

    % ---- (b) f_ad ----
    ax = nexttile(tl); hold(ax, 'on')
    for jj = 1:numel(idx)
        plot(ax, [lowest_base.f_ad_used(jj), lowest_base.f_ad_intercept(jj)], [y_pos(jj), y_pos(jj)], '-', ...
            'Color', style.gray_light, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
    plot(ax, lowest_base.f_ad_used, y_pos, 's', 'MarkerSize', 8, ...
        'MarkerFaceColor', style.lowest, 'MarkerEdgeColor', style.lowest);
    plot(ax, lowest_base.f_ad_intercept, y_pos, '^', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', style.fit, 'LineWidth', 1.5);
    xline(ax, 1, 'Color', style.gray_mid, 'LineWidth', 1)
    yticks(ax, 1:numel(idx)); yticklabels(ax, {})
    ylim(ax, [0.4, numel(idx) + 0.6])
    xlabel(ax, 'f_{ad} = LWP_{obs} / LWP_{ad}')
    ratio = lowest_base.f_ad_used ./ lowest_base.f_ad_intercept;
    title(ax, sprintf('(b) f_{ad} used / f_{ad} with the intercept as base: median %.2f, range %.2f-%.2f', ...
        median(ratio), min(ratio), max(ratio)), 'FontWeight', 'normal')
    grid(ax, 'on'); ax.GridAlpha = 0.1; ax.YGrid = 'off';

    title(tl, 'Profiles whose fitted line reaches LWC = 0 above the lowest kept sample: that sample is the base')

    finish_figure(fig, inputs, 'wood_adiab_fig4_lowest_sample_base')


    %% Figure 5: the LWC profiles behind Figure 4

    n_cols = min(5, numel(idx));
    n_rows = ceil(numel(idx) / n_cols);

    fig = make_figure(inputs, 'Lowest sample as the base: profiles', [30, 30, 360 * n_cols, 330 * n_rows + 120]);
    tl = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for jj = 1:numel(idx)
        ax = nexttile(tl);
        draw_profile_panel(ax, profiles(idx(jj)), style, 'above');
        if jj == 1
            ax_legend = ax;
        end
    end

    add_panel_legend(ax_legend, style, false, true);
    xlabel(tl, 'LWC (g m^{-3})')
    ylabel(tl, 'height above the lowest kept sample (m)')
    title(tl, ['The line reaches LWC = 0 above the lowest sample, so the lowest sample is the base. ', ...
        'Open circles: samples below where the line reaches zero, which count toward LWP_{obs}'])

    finish_figure(fig, inputs, 'wood_adiab_fig5_lowest_sample_base_profiles')

end


%% Figure 6: mean in-cloud T and p, and the Gamma_ad computed from them

fig = make_figure(inputs, 'Mean in-cloud T, p and Gamma_ad', [40, 40, 1600, 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
marker = {'o', '^'};

% ---- (a, b, c) distributions of T, p and Gamma_ad ----------------------
quantities = {summary.T_mean_C, summary.P_mean_hPa, summary.Gamma_ad_g_m3_per_km};
x_labels = {'mean in-cloud temperature (\circC)', 'mean in-cloud pressure (hPa)', ...
    '\Gamma_{ad} (g m^{-3} km^{-1})'};
panel_titles = {'(a) Temperature', '(b) Pressure', '(c) \Gamma_{ad}'};
bin_width = [1, 10, 0.05];
number_format = {'%.1f', '%.0f', '%.2f'};          % T in deg C, p in hPa, Gamma_ad in g/m^3/km

for qq = 1:3

    ax = nexttile(tl); hold(ax, 'on')
    values = quantities{qq};
    edges = floor(min(values) / bin_width(qq)) * bin_width(qq) : bin_width(qq) : ...
        ceil(max(values) / bin_width(qq)) * bin_width(qq) + bin_width(qq);
    y_max = 0;
    h = gobjects(0); names = {};

    if qq == 3
        % Reference values for marine stratus and stratocumulus: usually 2.0-2.8, and
        % 2.1-2.6 with a median of 2.24 over the southeast Pacific in EPIC 2001
        % (Zhou et al. 2006, GRL 33, L02807, doi:10.1029/2005GL024796, para. 8, citing
        % Albrecht et al. 1990, Austin et al. 1995 and Wood and Taylor 2001)
        h(end+1) = patch(ax, [2.0, 2.8, 2.8, 2.0], [0, 0, 1, 1] .* 1e3, style.gray_mid, ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none'); %#ok<SAGROW>
        names{end+1} = 'usual range 2.0-2.8 (Zhou et al. 2006)'; %#ok<SAGROW>
        h(end+1) = xline(ax, 2.24, 'Color', style.gray_dark, 'LineWidth', 1); %#ok<SAGROW>
        names{end+1} = 'EPIC 2001 median 2.24 (range 2.1-2.6)'; %#ok<SAGROW>
    end

    for cc = 1:2
        k = summary.campaign == campaign_names(cc) & isfinite(values);
        if ~any(k), continue, end
        counts = histcounts(values(k), edges);
        h(end+1) = stairs(ax, edges, [counts, counts(end)], 'Color', style.campaign(cc, :), 'LineWidth', 2); %#ok<SAGROW>
        fmt = number_format{qq};
        names{end+1} = sprintf(['%s: median ', fmt, ', range ', fmt, ' to ', fmt], campaign_names(cc), ...
            median(values(k)), min(values(k)), max(values(k))); %#ok<SAGROW>
        y_max = max(y_max, max(counts));
    end

    legend(ax, h, names, 'Location', 'northoutside', 'Box', 'off', 'FontSize', 8)
    xlim(ax, edges([1, end])); ylim(ax, [0, 1.1 * y_max])
    xlabel(ax, x_labels{qq}); ylabel(ax, 'profiles')
    title(ax, panel_titles{qq}, 'FontWeight', 'normal')
    grid(ax, 'on'); ax.GridAlpha = 0.1;

end

% ---- (d) Gamma_ad over the temperature-pressure plane -----------------
% The contours come from compute_adiabaticity_Wood2005 itself (gamma_ad_at below),
% so they are exactly the pipeline's thermodynamics, not a copy of it.
ax = nexttile(tl, [1, 2]); hold(ax, 'on')

T_grid_C = linspace(floor(min(summary.T_mean_C)) - 2, ceil(max(summary.T_mean_C)) + 2, 45);
P_grid_hPa = linspace(floor(min(summary.P_mean_hPa) / 10) * 10 - 20, ceil(max(summary.P_mean_hPa) / 10) * 10 + 20, 45);
[T_mesh, P_mesh] = meshgrid(T_grid_C, P_grid_hPa);
G_mesh = arrayfun(@gamma_ad_at, T_mesh, P_mesh);

[C, hc] = contour(ax, T_mesh, P_mesh, G_mesh, 1.4:0.1:3.0, 'LineColor', style.gray_mid, 'LineWidth', 0.75);
clabel(C, hc, 'FontSize', 8, 'Color', style.gray_dark, 'LabelSpacing', 300)

h = gobjects(0); names = {};
for cc = 1:2
    k = summary.campaign == campaign_names(cc) & has_gamma;
    if ~any(k), continue, end
    h(end+1) = scatter(ax, summary.T_mean_C(k), summary.P_mean_hPa(k), 30, style.campaign(cc, :), ...
        marker{cc}, 'filled', 'MarkerFaceAlpha', 0.8); %#ok<SAGROW>
    names{end+1} = sprintf('%s (%d profiles)', campaign_names(cc), sum(k)); %#ok<SAGROW>
end

set(ax, 'YDir', 'reverse')
xlim(ax, T_grid_C([1, end])); ylim(ax, P_grid_hPa([1, end]))
legend(ax, h, names, 'Location', 'southwest', 'Box', 'off')
xlabel(ax, 'mean in-cloud temperature (\circC)')
ylabel(ax, 'mean in-cloud pressure (hPa)')
title(ax, '(d) Each profile''s mean in-cloud state; gray contours: \Gamma_{ad} (g m^{-3} km^{-1})', ...
    'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1;

% ---- (e) the equations ------------------------------------------------
ax = nexttile(tl); axis(ax, 'off')
equations = { ...
    '\textbf{How $\Gamma_{ad}$ is computed} (compute\_adiabaticity\_Wood2005)', ...
    '', ...
    '$\displaystyle \Gamma_{ad} = \rho_{air}\,\frac{c_p}{L_v}\,\left(\Gamma_d - \Gamma_m\right)$', ...
    '', ...
    '$\displaystyle \rho_{air} = \frac{p}{R_d\,T}, \qquad \Gamma_d = \frac{g}{c_p}$', ...
    '', ...
    '$\displaystyle \Gamma_m = g\,\frac{1 + \frac{L_v\,q_s}{R_d\,T}}{c_p + \frac{L_v^2\,q_s\,\epsilon}{R_d\,T^2}}$ \quad (Randall 2009, eq. 36)', ...
    '', ...
    '$\displaystyle q_s = \frac{\epsilon\,e_s}{p - (1 - \epsilon)\,e_s}, \qquad e_s = 611.2\,\exp\left(\frac{17.67\,T_C}{T_C + 243.5}\right)$ Pa', ...
    '\quad (Bolton 1980)', ...
    '', ...
    '$g = 9.81$ m s$^{-2}$, \ $c_p = 1005$, \ $R_d = 287.04$, \ $R_v = 461.5$ J kg$^{-1}$ K$^{-1}$', ...
    '$L_v = 2.50 \times 10^{6}$ J kg$^{-1}$, \ $\epsilon = R_d / R_v$', ...
    '', ...
    '$T$, $p$: trapezoid-weighted mean over the sampled cloud depth', ...
    '(Wood 2005: one $\Gamma_{ad}$ per cloud, at the mean in-cloud state).', ...
    'Units: kg m$^{-3}$ per m; $\times 10^{6}$ gives g m$^{-3}$ km$^{-1}$.'};
text(ax, 0, 1, equations, 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontSize', 12)

title(tl, 'One \Gamma_{ad} per profile, at its mean in-cloud temperature and pressure')

finish_figure(fig, inputs, 'wood_adiab_fig6_gamma_ad')


%% Figure 7: climbs against descents
% Whether the aircraft flew up or down through a cloud should not change the
% cloud's adiabaticity. A difference has to come from (i) the in-cloud state used
% for Gamma_ad, (ii) the cloud base and thickness, (iii) the LWC measurement itself,
% or (iv) climbs and descents sampling different clouds. Each panel compares one
% ingredient. Bootstrap 95% CIs on the climb-minus-descent median difference.

rand_stream = RandStream('mt19937ar', 'Seed', inputs.rng_seed);
fprintf('\nClimb-descent bootstrap: seed %d, %d resamples\n', inputs.rng_seed, inputs.n_bootstrap);

quantities = {summary.T_mean_C, summary.P_mean_hPa, summary.Gamma_ad_g_m3_per_km, ...
    summary.f_ad, summary.h_m, -summary.z_cb_above_lowest_m, ...
    summary.lwp_obs_g_m2, summary.abs_w_m_s, summary.slope_over_gamma_ad};
y_labels = {'mean in-cloud T (\circC)', 'mean in-cloud p (hPa)', '\Gamma_{ad} (g m^{-3} km^{-1})', ...
    'f_{ad}', 'cloud thickness h (m)', 'z_{cb} below the lowest sample (m)', ...
    'LWP_{obs} (g m^{-2})', '|vertical speed| (m s^{-1})', 'fitted slope / \Gamma_{ad}'};
diff_format = {'%+.2f', '%+.1f', '%+.3f', '%+.2f', '%+.0f', '%+.1f', '%+.0f', '%+.2f', '%+.2f'};
panel_titles = {'(a) Temperature', '(b) Pressure', '(c) \Gamma_{ad}', '(d) f_{ad}', '(e) Thickness', ...
    '(f) Extrapolation below the lowest sample', '(g) LWP_{obs}', '(h) Climb or descent rate', ...
    '(i) Bottom-third slope'};

fig = make_figure(inputs, 'Climbs vs descents', [40, 40, 1500, 1050]);
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

direction_table = table();

for qq = 1:numel(quantities)

    ax = nexttile(tl); hold(ax, 'on')
    values = quantities{qq};
    tick_pos = []; tick_lab = {};
    diff_text = {};

    for cc = 1:2
        for dd = 1:2                                       % 1 = climb, 2 = descent
            k = summary.campaign == campaign_names(cc) & has_gamma & ...
                summary.is_ascending == (dd == 1) & isfinite(values);
            if ~any(k), continue, end
            x0 = 3 * (cc - 1) + dd;
            jitter = 0.28 .* (rand(rand_stream, sum(k), 1) - 0.5);
            if dd == 1
                scatter(ax, x0 + jitter, values(k), 14, style.campaign(cc, :), 'filled', 'MarkerFaceAlpha', 0.6)
            else
                scatter(ax, x0 + jitter, values(k), 14, style.campaign(cc, :), 'LineWidth', 0.8)
            end
            plot(ax, x0 + [-0.25, 0.25], median(values(k)) .* [1, 1], '-', 'Color', style.gray_dark, 'LineWidth', 2)
            tick_pos(end+1) = x0; %#ok<SAGROW>
            tick_lab{end+1} = sprintf('%s\\newline%s', campaign_short(cc), direction_word(dd == 1)); %#ok<SAGROW>
        end

        k_up   = summary.campaign == campaign_names(cc) & has_gamma & summary.is_ascending;
        k_down = summary.campaign == campaign_names(cc) & has_gamma & ~summary.is_ascending;
        if any(k_up) && any(k_down)
            [d, ci] = median_difference_ci(values(k_up), values(k_down), rand_stream, inputs.n_bootstrap);
            fmt = diff_format{qq};
            diff_text{end+1} = sprintf(['%s: \\Delta = ', fmt, ' [', fmt, ', ', fmt, ']'], ...
                campaign_short(cc), d, ci(1), ci(2)); %#ok<SAGROW>
            direction_table = [direction_table; table(string(y_labels{qq}), campaign_names(cc), ...
                median(values(k_up), 'omitnan'), median(values(k_down), 'omitnan'), d, ci(1), ci(2), ...
                'VariableNames', {'quantity', 'campaign', 'median_climb', 'median_descent', ...
                'climb_minus_descent', 'ci_low', 'ci_high'})]; %#ok<AGROW>
        end
    end

    xticks(ax, tick_pos); xticklabels(ax, tick_lab)
    xlim(ax, [0.4, 5.6])

    % A few extreme profiles (the 30-Aug-2017 ORACLES pair, a drizzle shaft) would
    % squash everything else, so the axis spans the central 99% and says what it hides
    plotted = values(has_gamma & isfinite(values));
    y_range = prctile(plotted, [0.5, 99.5]);
    y_pad = 0.08 .* diff(y_range);
    ylim(ax, y_range + [-y_pad, y_pad])
    n_hidden = sum(plotted < y_range(1) - y_pad | plotted > y_range(2) + y_pad);
    if n_hidden > 0
        text(ax, 0.98, 0.97, sprintf('%d off scale', n_hidden), 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', style.gray_dark)
    end
    ylabel(ax, y_labels{qq})
    title(ax, {panel_titles{qq}, strjoin(diff_text, '   ')}, 'FontWeight', 'normal', 'FontSize', 9)
    grid(ax, 'on'); ax.GridAlpha = 0.1; ax.XGrid = 'off';

end

title(tl, sprintf(['Climbs (filled) vs descents (open); bar = median; \\Delta = climb - descent median ', ...
    'with bootstrap 95%% CI (seed %d)'], inputs.rng_seed))

disp([newline, 'Climb minus descent, median and bootstrap 95% CI:'])
disp(direction_table)

finish_figure(fig, inputs, 'wood_adiab_fig7_climb_vs_descent')


%% Figure 8: where in the cloud do climbs and descents differ?
% The per-level ratio LWC(z) / q_ad(z) against z* = (z - z_cb) / h (Wood's Fig. 4c),
% median in z* bins. A difference only near z* = 0 points at the cloud base; one
% near the top at entrainment or cloud-top sampling; one at every level at the
% LWC measurement or the thickness.

fig = make_figure(inputs, 'Climbs vs descents: per-level ratio', [40, 40, 1100, 560]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
z_edges = 0:0.1:1;
z_centers = z_edges(1:end-1) + 0.05;

for cc = 1:2

    ax = nexttile(tl); hold(ax, 'on')
    h = gobjects(0); names = {};

    for dd = 1:2
        k = find(summary.campaign == campaign_names(cc) & has_gamma & summary.is_ascending == (dd == 1));
        if isempty(k), continue, end
        z_all = cell2mat(arrayfun(@(p) p.a.z_star, profiles(k), 'UniformOutput', false));
        f_all = cell2mat(arrayfun(@(p) p.a.f_ad_z, profiles(k), 'UniformOutput', false));
        f_med = nan(size(z_centers)); f_q1 = f_med; f_q3 = f_med;
        for bb = 1:numel(z_centers)
            in_bin = z_all >= z_edges(bb) & z_all < z_edges(bb + 1) & isfinite(f_all);
            if sum(in_bin) >= 5
                f_med(bb) = median(f_all(in_bin)); f_q1(bb) = prctile(f_all(in_bin), 25); f_q3(bb) = prctile(f_all(in_bin), 75);
            end
        end
        ok = isfinite(f_q1);
        patch(ax, [f_q1(ok), fliplr(f_q3(ok))], [z_centers(ok), fliplr(z_centers(ok))], style.campaign(cc, :), ...
            'FaceAlpha', 0.12 + 0.08 * (dd == 1), 'EdgeColor', 'none', 'HandleVisibility', 'off')
        line_style = {'-', '--'};
        h(end+1) = plot(ax, f_med, z_centers, line_style{dd}, 'Color', style.campaign(cc, :), 'LineWidth', 2); %#ok<SAGROW>
        names{end+1} = sprintf('%ss (%d), median and 25-75%%', direction_word(dd == 1), numel(k)); %#ok<SAGROW>
    end

    xline(ax, 1, 'Color', style.gray_mid, 'LineWidth', 1, 'HandleVisibility', 'off')
    xlim(ax, [0, 2]); ylim(ax, [0, 1])
    legend(ax, h, names, 'Location', 'northeast', 'Box', 'off')
    xlabel(ax, 'LWC(z) / q_{ad}(z)')
    ylabel(ax, 'z_* = (z - z_{cb}) / h')
    title(ax, campaign_names(cc), 'FontWeight', 'normal')
    grid(ax, 'on'); ax.GridAlpha = 0.1;

end

title(tl, 'Per-level adiabaticity by flight direction (line-intercept base)')
finish_figure(fig, inputs, 'wood_adiab_fig8_climb_vs_descent_levels')


%% Figure 9: two tests of where the climb-descent difference comes from
% (a) NEIGHBORING PAIRS. A climb and a descent flown within
%     inputs.pair_max_separation_s of each other on the same flight sample nearly the
%     same cloud. If the difference survives inside pairs it is a measurement or
%     processing effect; if it vanishes, climbs and descents sampled different clouds.
% (b) TIME SHIFT BETWEEN LWC AND ALTITUDE. If the LWC record lagged the altitude record
%     by tau seconds, every LWC value would be plotted tau*w higher than it was measured
%     on a climb and tau*w lower on a descent. That makes climbs look thinner, with the
%     base too high, and descents thicker, which is the observed sign. Each profile's
%     f_ad is recomputed with the LWC record shifted by tau seconds (positive tau pairs
%     each altitude with the LWC recorded tau s later, which undoes a lag of tau s).
%     A zero crossing near a small positive tau is consistent with a lag; it does not
%     prove one, because any direction-dependent bias can be cancelled by some shift.
%     ORACLES is masked outside cloud, so a shift also pulls zeros into one end of the
%     profile, which exaggerates its sensitivity.
% (c) CLOUD TOP IN PAIRS. Cloud top is sharp, so it is the cleanest place to look for
%     a shift. A lag of tau puts the top 2 w tau higher on the climb than on the
%     neighboring descent. Something that lingers after the aircraft leaves the cloud
%     (water on a probe, for example) would do the same, so this test cannot tell the
%     two apart. On VOCALS-REx the finder's start/stop asymmetry also lets climbs run
%     a little further past the top; on ORACLES the mask removes that effect.

fig = make_figure(inputs, 'Climb-descent tests', [40, 40, 1750, 560]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) neighboring climb-descent pairs ----
ax = nexttile(tl); hold(ax, 'on')
pair_text = {};

for cc = 1:2
    [up, down] = climb_descent_pairs(summary, campaign_names(cc), inputs.pair_max_separation_s);
    paired_diff = summary.f_ad(up) - summary.f_ad(down);
    if isempty(paired_diff), continue, end
    jitter = 0.3 .* (rand(rand_stream, numel(paired_diff), 1) - 0.5);
    scatter(ax, cc + jitter, paired_diff, 26, style.campaign(cc, :), 'filled', 'MarkerFaceAlpha', 0.7)
    boot = median(paired_diff(randi(rand_stream, numel(paired_diff), numel(paired_diff), inputs.n_bootstrap)), 1);
    ci = prctile(boot, [2.5, 97.5]);
    plot(ax, cc + [-0.25, 0.25], median(paired_diff) .* [1, 1], '-', 'Color', style.gray_dark, 'LineWidth', 2)
    pair_text{end+1} = sprintf('%s: %d pairs, median %+.2f [%+.2f, %+.2f]', campaign_names(cc), ...
        numel(paired_diff), median(paired_diff), ci(1), ci(2)); %#ok<SAGROW>
end

yline(ax, 0, 'Color', style.gray_mid, 'LineWidth', 1)
xticks(ax, 1:2); xticklabels(ax, campaign_names); xlim(ax, [0.4, 2.6])
ylabel(ax, 'f_{ad}(climb) - f_{ad}(descent) within a pair')
title(ax, [{sprintf('(a) Climb-descent pairs less than %.0f min apart', inputs.pair_max_separation_s / 60)}, pair_text], ...
    'FontWeight', 'normal', 'FontSize', 9)
grid(ax, 'on'); ax.GridAlpha = 0.1; ax.XGrid = 'off';
disp([newline, strjoin(pair_text, newline)])

% ---- (b) time shift between the LWC and altitude records ----
ax = nexttile(tl); hold(ax, 'on')
shifts = inputs.lwc_time_shifts_s;
h = gobjects(0); names = {};

for cc = 1:2
    k_up   = find(summary.campaign == campaign_names(cc) & has_gamma & summary.is_ascending);
    k_down = find(summary.campaign == campaign_names(cc) & has_gamma & ~summary.is_ascending);
    if isempty(k_up) || isempty(k_down), continue, end
    F_up = vertcat(profiles(k_up).f_ad_by_shift);       % [n_climbs x n_shifts]
    F_down = vertcat(profiles(k_down).f_ad_by_shift);   % [n_descents x n_shifts]
    d = nan(size(shifts)); ci = nan(2, numel(shifts));
    for ss = 1:numel(shifts)
        [d(ss), ci(:, ss)] = median_difference_ci(F_up(:, ss), F_down(:, ss), rand_stream, inputs.n_bootstrap);
    end
    patch(ax, [shifts, fliplr(shifts)], [ci(1, :), fliplr(ci(2, :))], style.campaign(cc, :), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    h(end+1) = plot(ax, shifts, d, '-o', 'Color', style.campaign(cc, :), 'LineWidth', 2, ...
        'MarkerFaceColor', style.campaign(cc, :), 'MarkerSize', 5); %#ok<SAGROW>
    names{end+1} = char(campaign_names(cc)); %#ok<SAGROW>
    fprintf('%s time-shift scan, climb - descent median f_ad:', campaign_names(cc));
    fprintf(' %+d s: %+.3f |', [shifts; d]); fprintf('\n');
end

yline(ax, 0, 'Color', style.gray_mid, 'LineWidth', 1)
xline(ax, 0, 'Color', style.gray_mid, 'LineWidth', 0.5)
legend(ax, h, names, 'Location', 'northeast', 'Box', 'off')
xlabel(ax, 'LWC record shifted by \tau (s); \tau > 0 undoes an LWC lag of \tau seconds')
ylabel(ax, 'median f_{ad}: climbs - descents')
title(ax, '(b) Climb-descent difference if the LWC record were shifted in time (95% CI shaded)', ...
    'FontWeight', 'normal', 'FontSize', 9)
grid(ax, 'on'); ax.GridAlpha = 0.1;

% ---- (c) cloud top and lowest sample in neighboring pairs ----
ax = nexttile(tl); hold(ax, 'on')
top_text = {};

for cc = 1:2
    [up, down] = climb_descent_pairs(summary, campaign_names(cc), inputs.pair_max_separation_s);
    if isempty(up), continue, end
    d_top = summary.z_top_m(up) - summary.z_top_m(down);              % m
    d_low = summary.z_lowest_m(up) - summary.z_lowest_m(down);        % m
    w_pair = 0.5 .* (summary.abs_w_m_s(up) + summary.abs_w_m_s(down)); % m/s
    for jj = 1:2
        values = d_top; x0 = cc; if jj == 2, values = d_low; x0 = cc + 3; end
        jitter = 0.3 .* (rand(rand_stream, numel(values), 1) - 0.5);
        scatter(ax, x0 + jitter, values, 22, style.campaign(cc, :), 'filled', 'MarkerFaceAlpha', 0.6)
        plot(ax, x0 + [-0.25, 0.25], median(values) .* [1, 1], '-', 'Color', style.gray_dark, 'LineWidth', 2)
    end
    boot = median(d_top(randi(rand_stream, numel(d_top), numel(d_top), inputs.n_bootstrap)), 1);
    ci = prctile(boot, [2.5, 97.5]);
    top_text{end+1} = sprintf('%s top: %+.1f m [%+.1f, %+.1f], as a lag %.1f s', campaign_short(cc), ...
        median(d_top), ci(1), ci(2), median(d_top) ./ (2 .* median(w_pair))); %#ok<SAGROW>
end

yline(ax, 0, 'Color', style.gray_mid, 'LineWidth', 1)
xticks(ax, [1, 2, 4, 5]); xticklabels(ax, {'VOCALS\newlinecloud top', 'ORACLES\newlinecloud top', ...
    'VOCALS\newlinelowest sample', 'ORACLES\newlinelowest sample'})
xlim(ax, [0.4, 5.6]); ylim(ax, [-150, 150])
ylabel(ax, 'climb - descent within a pair (m)')
title(ax, [{'(c) Cloud top and lowest-sample height, climb - descent'}, top_text], 'FontWeight', 'normal', 'FontSize', 9)
grid(ax, 'on'); ax.GridAlpha = 0.1; ax.XGrid = 'off';
disp(strjoin(top_text, newline))

finish_figure(fig, inputs, 'wood_adiab_fig9_climb_descent_tests')




%% ========================================================================
%  Local functions
%  ========================================================================

function r = profile_record(campaign_name, label, flight_id, vert_prof)
% The pieces of one profile that the figures need

r.campaign = campaign_name;
r.label = label;
r.flight = flight_id;
r.t_start_s = double(vert_prof.time(1));                                   % s
r.w_m_s = mean(diff(double(vert_prof.altitude)) ./ diff(double(vert_prof.time)));   % m/s, mean vertical speed
r.is_ascending = r.w_m_s > 0;
r.z_m = double(vert_prof.altitude(:))';                     % m, sampling order
r.lwc_g_m3 = double(vert_prof.lwc(:))';                     % g/m^3, total LWC as used for f_ad
r.a = vert_prof.adiabaticity;                               % default: line-intercept base
r.a_lowest = vert_prof.adiabaticity_lowest_sample;          % lowest sample taken as the base
r.f_ad_by_shift = [];                                       % filled by f_ad_with_time_shift

end



function f_ad = f_ad_with_time_shift(time_s, altitude_m, lwc_g_m3, temp_C, pres_hPa, ...
    profile_time_s, lwc_threshold_g_m3, shifts_s)
% f_ad (line-intercept base) of one profile after pairing each altitude with the LWC
% recorded shift_s seconds later, for every shift in shifts_s. A positive shift undoes
% an LWC record that lags the altitude record by that many seconds. Temperature and
% pressure stay with the altitude. The profile's own samples are found in the full
% flight time series by matching its first and last time stamps.

i_start = find(time_s == profile_time_s(1), 1);
i_end   = find(time_s == profile_time_s(end), 1);
idx = i_start:i_end;
f_ad = nan(1, numel(shifts_s));

for kk = 1:numel(shifts_s)

    idx_lwc = idx + shifts_s(kk);
    if idx_lwc(1) < 1 || idx_lwc(end) > numel(lwc_g_m3)
        continue
    end

    a = compute_adiabaticity_Wood2005(altitude_m(idx), lwc_g_m3(idx_lwc), temp_C(idx), pres_hPa(idx), ...
        'base_method', 'extrapolate', 'lwc_threshold_g_m3', lwc_threshold_g_m3);
    f_ad(kk) = a.f_ad_lwp;

end

end



function [f_ad, h_m, lwp_obs_g_m2, lwp_ad_g_m2] = f_ad_with_intercept_base(r)
% f_ad if the base were put where the fitted line reaches LWC = 0, even though that
% is above the lowest sample. The samples below that height are left out of LWP_obs.
% Same line and Gamma_ad as compute_adiabaticity_Wood2005; only for profiles whose
% base_source is 'lowest sample, intercept above it'.

a = r.a;
Gamma_g_m3_per_m = a.Gamma_ad_g_m3_per_km ./ 1e3;                          % g/m^3 per m
z_line_zero_m = a.z_lowest_sample_m - a.lwc_fit_at_lowest_g_m3 ./ (a.fit_slope_g_m3_per_km ./ 1e3);   % m

above = isfinite(r.lwc_g_m3) & r.z_m >= z_line_zero_m;
h_m = a.z_top_m - z_line_zero_m;                                            % m
lwp_obs_g_m2 = abs(trapz(r.z_m(above), r.lwc_g_m3(above)));                 % g/m^2
lwp_ad_g_m2 = 0.5 .* Gamma_g_m3_per_m .* h_m.^2;                            % g/m^2
f_ad = lwp_obs_g_m2 ./ lwp_ad_g_m2;

end



function [up, down] = climb_descent_pairs(summary, campaign_name, max_separation_s)
% Neighboring climb-descent pairs: consecutive profiles on the same flight, flown in
% opposite directions, starting less than max_separation_s apart. Each profile is used
% at most once. Returns row indices into summary for the climb and the descent.

up = []; down = [];
flights = unique(summary.flight(summary.campaign == campaign_name));

for ff = 1:numel(flights)

    k = find(summary.campaign == campaign_name & summary.flight == flights(ff) & isfinite(summary.f_ad));
    [~, order] = sort(summary.t_start_s(k));
    k = k(order);

    jj = 1;
    while jj < numel(k)
        pair = k([jj, jj + 1]);
        if summary.is_ascending(pair(1)) ~= summary.is_ascending(pair(2)) && ...
                summary.t_start_s(pair(2)) - summary.t_start_s(pair(1)) < max_separation_s
            up(end+1, 1)   = pair(summary.is_ascending(pair));  %#ok<AGROW>
            down(end+1, 1) = pair(~summary.is_ascending(pair)); %#ok<AGROW>
            jj = jj + 2;
        else
            jj = jj + 1;
        end
    end

end

end



function [d, ci] = median_difference_ci(x, y, rand_stream, n_bootstrap)
% Difference of medians, median(x) - median(y), with a bootstrap 95% confidence
% interval. x and y are resampled independently, with replacement.

x = x(isfinite(x)); x = x(:);
y = y(isfinite(y)); y = y(:);
d = median(x) - median(y);

boot_x = median(x(randi(rand_stream, numel(x), numel(x), n_bootstrap)), 1);   % [1 x n_bootstrap]
boot_y = median(y(randi(rand_stream, numel(y), numel(y), n_bootstrap)), 1);   % [1 x n_bootstrap]
ci = prctile(boot_x - boot_y, [2.5, 97.5]);

end



function reason = reason_for_fallback(a)
% Why compute_adiabaticity_Wood2005 rejected the line fit, in the same order it tests

reason = "";

if ~strcmp(a.base_source, 'fallback')
    return
end

slope_ratio = a.fit_slope_over_gamma_ad;
limits = a.slope_limits_over_gamma_ad;

if a.n_fit_samples < 3 || ~isfinite(slope_ratio)
    reason = "fewer than 3 samples to fit";
elseif slope_ratio <= 0
    reason = sprintf("slope %.2f Gamma_ad <= 0", slope_ratio);
elseif slope_ratio <= limits(1)
    reason = sprintf("slope %.2f Gamma_ad <= %.2f", slope_ratio, limits(1));
else
    reason = sprintf("slope %.2f Gamma_ad > %g", slope_ratio, limits(2));
end

end



function word = direction_word(is_ascending)
% "climb" or "descent" for each profile (string array, same size as the input)

word = repmat("descent", size(is_ascending));
word(is_ascending) = "climb";

end



function draw_profile_panel(ax, r, style, mode)
% One profile: LWC against height above its lowest kept sample, with the fitted line,
% the cloud base it gives, and the adiabatic profile from that base.
%
% mode: 'fit' (every profile), 'fallback' (adds the reason for the fallback) or
%       'above' (adds the lowest sample as an alternative base)

a = r.a;
z_m = r.z_m;
lwc = r.lwc_g_m3;

hold(ax, 'on')
ax.FontSize = 7;
grid(ax, 'on'); ax.GridAlpha = 0.1;

if r.is_ascending
    direction = 'climb';
else
    direction = 'descent';
end

% Two ORACLES profiles have no temperature, so no Gamma_ad and no base
if ~isfinite(a.Gamma_ad_g_m3_per_km)
    plot(ax, lwc, z_m - min(z_m), '.', 'Color', style.gray_light, 'MarkerSize', 6)
    title(ax, sprintf('%s (%s)', r.label, direction), 'FontWeight', 'normal', 'FontSize', 8)
    text(ax, 0.05, 0.9, 'no temperature: f_{ad} not computed', 'Units', 'normalized', 'FontSize', 7)
    return
end

z_lowest_m = a.z_lowest_sample_m;
h_sampled_m = a.h_lowest_sample_m;
y_m = z_m - z_lowest_m;                                    % m above the lowest sample

% The same samples compute_adiabaticity_Wood2005 fits
in_fit = isfinite(lwc) & y_m ./ h_sampled_m <= a.fit_fraction + 1e-4;

Gamma_g_m3_per_m = a.Gamma_ad_g_m3_per_km ./ 1e3;           % g/m^3 per m
slope_g_m3_per_m = a.fit_slope_g_m3_per_km ./ 1e3;          % g/m^3 per m
lwc_line_at_lowest = a.lwc_fit_at_lowest_g_m3;              % g/m^3
y_cb_m = 0 - a.dz_base_m;                                   % base, m above the lowest sample (<= 0; 0 - keeps a zero from printing as -0)
y_top_fit_m = a.fit_fraction .* h_sampled_m;                % top of the fitted layer
is_fallback = strcmp(a.base_source, 'fallback');
is_lowest_base = strcmp(a.base_source, 'lowest sample, intercept above it');

% Where the fitted line reaches LWC = 0, m above the lowest sample (positive only
% when the line reaches zero above it and the lowest sample is the base)
y_line_zero_m = -lwc_line_at_lowest ./ slope_g_m3_per_m;

% Axis limits: from just below the base (or the lowest sample) to the cloud top
y_low = min(0, y_cb_m);
y_low = y_low - 0.04 .* (h_sampled_m - y_low);
y_high = 1.02 .* h_sampled_m;
x_high = 1.08 .* max([lwc, Gamma_g_m3_per_m .* (h_sampled_m - y_cb_m), lwc_line_at_lowest]);

% Samples: above the fitted layer in gray, the fitted ones in the fit color
plot(ax, lwc(~in_fit), y_m(~in_fit), '.', 'Color', style.gray_light, 'MarkerSize', 6)
plot(ax, lwc(in_fit), y_m(in_fit), '.', 'Color', style.fit, 'MarkerSize', 8)
yline(ax, y_top_fit_m, 'Color', style.gray_mid, 'LineWidth', 0.5)

% Adiabatic LWC from the base that was used: Gamma_ad (z - z_cb)
plot(ax, Gamma_g_m3_per_m .* ([y_cb_m, h_sampled_m] - y_cb_m), [y_cb_m, h_sampled_m], '-', ...
    'Color', style.gray_dark, 'LineWidth', 1)

if isfinite(slope_g_m3_per_m)

    if ~is_fallback

        % The fitted line, from where it reaches LWC = 0 (or the lowest sample) to the
        % top of the fitted layer. Solid over the data it was fit to, dashed where it
        % is extrapolated below the lowest sample.
        y_line = [max(0, y_cb_m), y_top_fit_m];
        if is_lowest_base
            y_line(1) = y_line_zero_m;
        end
        plot(ax, lwc_line_at_lowest + slope_g_m3_per_m .* y_line, y_line, '-', 'Color', style.fit, 'LineWidth', 1.5)
        if y_cb_m < 0
            plot(ax, [0, lwc_line_at_lowest], [y_cb_m, 0], '--', 'Color', style.fit, 'LineWidth', 1.5)
        end

    else

        % The rejected line, carried down to the bottom of the axis to show where it heads
        y_line = [y_low, y_top_fit_m];
        x_line = lwc_line_at_lowest + slope_g_m3_per_m .* y_line;
        plot(ax, x_line, y_line, ':', 'Color', style.fit, 'LineWidth', 1.5)

    end

end

if is_fallback
    % Fallback: an adiabat from the LWC threshold at the lowest sample down to LWC = 0
    plot(ax, [0, a.lwc_threshold_g_m3], [y_cb_m, 0], '-', 'Color', style.fallback, 'LineWidth', 2)
    plot(ax, 0, y_cb_m, '^', 'MarkerSize', 7, 'MarkerFaceColor', style.fallback, 'MarkerEdgeColor', 'w')
elseif is_lowest_base
    % The line reaches LWC = 0 above the lowest sample (open triangle, not used);
    % the lowest sample is the base (square)
    plot(ax, 0, y_line_zero_m, '^', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', style.fit, 'LineWidth', 1.2)
    plot(ax, 0, 0, 's', 'MarkerSize', 7, 'MarkerFaceColor', style.lowest, 'MarkerEdgeColor', 'w')
    if strcmp(mode, 'above')
        % Samples between the lowest sample and where the line reaches zero
        below = y_m < y_line_zero_m;
        plot(ax, lwc(below), y_m(below), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', style.lowest, 'LineWidth', 1)
    end
else
    plot(ax, 0, y_cb_m, '^', 'MarkerSize', 7, 'MarkerFaceColor', style.fit, 'MarkerEdgeColor', 'w')
end

xlim(ax, [0, x_high]); ylim(ax, [y_low, y_high])

% Panel text: the numbers that decide the base and the result
if is_fallback
    base_word = 'fallback';
elseif is_lowest_base
    base_word = 'lowest sample';
else
    base_word = 'intercept';
end
panel_text = sprintf('slope = %.2f \\Gamma_{ad}\nz_{cb} (%s): %+.0f m\nf_{ad} = %.2f', ...
    slope_g_m3_per_m ./ Gamma_g_m3_per_m, base_word, y_cb_m, a.f_ad_lwp);
if strcmp(mode, 'above')
    panel_text = sprintf('%s\nline reaches 0 at %+.1f m\nf_{ad}, intercept as base = %.2f', ...
        panel_text, y_line_zero_m, f_ad_with_intercept_base(r));
end
text(ax, 0.03, 0.97, panel_text, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7)

if strcmp(mode, 'fallback')
    reason = char(reason_for_fallback(a));
    if isfinite(slope_g_m3_per_m) && slope_g_m3_per_m > 0
        reason = sprintf('%s\n(line reaches 0 at %.0f m below)', reason, lwc_line_at_lowest ./ slope_g_m3_per_m);
    end
    reason = replace(reason, ["Gamma_ad", "<="], ["\Gamma_{ad}", "\leq"]);
    title(ax, {sprintf('%s (%s)', r.label, direction), char(reason)}, 'FontWeight', 'normal', 'FontSize', 8)
else
    title(ax, sprintf('%s (%s)', r.label, direction), 'FontWeight', 'normal', 'FontSize', 8)
end

end



function add_panel_legend(ax, style, include_fallback, include_lowest)
% One legend for a page of panels, built from invisible stand-ins so every entry
% appears whether or not the first panel uses it. ax is any axes in the page's
% tiled layout; the legend is moved to the layout's south tile.

h = gobjects(0); names = {};
h(end+1) = plot(ax, NaN, NaN, '.', 'Color', style.gray_light, 'MarkerSize', 12);     names{end+1} = 'sample above the fitted layer';
h(end+1) = plot(ax, NaN, NaN, '.', 'Color', style.fit, 'MarkerSize', 14);            names{end+1} = 'sample in the bottom third (fitted)';
h(end+1) = plot(ax, NaN, NaN, '-', 'Color', style.fit, 'LineWidth', 1.5);            names{end+1} = 'line fit';
h(end+1) = plot(ax, NaN, NaN, '--', 'Color', style.fit, 'LineWidth', 1.5);           names{end+1} = 'line extended to LWC = 0';
h(end+1) = plot(ax, NaN, NaN, '^', 'MarkerFaceColor', style.fit, 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
names{end+1} = 'cloud base z_{cb} (line intercept)';
h(end+1) = plot(ax, NaN, NaN, '-', 'Color', style.gray_dark, 'LineWidth', 1);        names{end+1} = 'adiabatic from z_{cb}: \Gamma_{ad} (z - z_{cb})';
h(end+1) = plot(ax, NaN, NaN, '-', 'Color', style.gray_mid, 'LineWidth', 0.5);       names{end+1} = 'top of the fitted layer';

if include_fallback
    h(end+1) = plot(ax, NaN, NaN, ':', 'Color', style.fit, 'LineWidth', 1.5);        names{end+1} = 'rejected line fit';
    h(end+1) = plot(ax, NaN, NaN, '-', 'Color', style.fallback, 'LineWidth', 2);     names{end+1} = 'fallback: adiabat from the LWC threshold';
    h(end+1) = plot(ax, NaN, NaN, '^', 'MarkerFaceColor', style.fallback, 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    names{end+1} = 'cloud base z_{cb} (fallback)';
end

if include_lowest
    h(end+1) = plot(ax, NaN, NaN, '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', style.fit, 'LineWidth', 1.2, 'MarkerSize', 8);
    names{end+1} = 'line reaches LWC = 0 above the lowest sample (not used)';
    h(end+1) = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', style.lowest, 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    names{end+1} = 'cloud base z_{cb} (lowest sample)';
    h(end+1) = plot(ax, NaN, NaN, 'o', 'MarkerEdgeColor', style.lowest, 'LineWidth', 1, 'MarkerSize', 6);
    names{end+1} = 'sample below where the line reaches 0';
end

lgd = legend(ax, h, names, 'NumColumns', 4, 'Box', 'off', 'FontSize', 9);
lgd.Layout.Tile = 'south';

end



function G = gamma_ad_at(T_C, P_hPa)
% Gamma_ad [g/m^3 per km] at one temperature and pressure, taken from
% compute_adiabaticity_Wood2005 itself by handing it a two-sample, isothermal,
% isobaric "profile". This keeps the contours identical to the pipeline's physics.

a = compute_adiabaticity_Wood2005([0, 1], [0, 0], [T_C, T_C], [P_hPa, P_hPa], 'base_method', 'lowest_sample');
G = a.Gamma_ad_g_m3_per_km;

end



function fig = make_figure(inputs, name, position)
% A figure with factory graphics defaults (new_factory_figure). Hidden when the
% figures are being saved, because exporting a visible figure can hang R2024b.

if isempty(inputs.save_figures_to)
    fig = new_factory_figure('Name', name, 'Position', position);
else
    fig = new_factory_figure('Name', name, 'Position', position, 'Visible', 'off');
end

end



function finish_figure(fig, inputs, file_stub)
% Save and close the figure if a folder was given; otherwise leave it on screen

if ~isempty(inputs.save_figures_to)
    if ~isfolder(inputs.save_figures_to)
        mkdir(inputs.save_figures_to)
    end
    exportgraphics(fig, fullfile(inputs.save_figures_to, [file_stub, '.png']), 'Resolution', 150)
    close(fig)
end

end
