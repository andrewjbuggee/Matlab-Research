%% Compare cloud-base definitions for the Wood (2005) adiabaticity
%
% DESCRIPTION:
%   f_ad = LWP_obs / LWP_ad depends on the cloud thickness squared, so how cloud
%   base is defined matters. There is no independent measurement of the true
%   base in these data, so this figure judges the definitions by evidence that
%   does not depend on them being right:
%
%   (a) Cumulative distribution of per-profile f_ad. A closed adiabatic parcel
%       cannot exceed 1, but a single slanted profile can, so this panel shows
%       the outcome and is not by itself a test of which definition is correct.
%
%   (b) How far below the lowest sample each profile's base is placed by the
%       measured definition and by the extrapolation (the default). Profiles
%       where the extrapolation used its fallback are marked.
%
%   (c) Median f_ad on climbs and descents, with bootstrap 95% confidence
%       intervals. Whether the aircraft flew up or down through a cloud should
%       not change the cloud's adiabaticity, so a definition that leaves a large
%       offset between the two has a sampling artifact in it. This is the most
%       direct test available here, but with so few descents it can only detect
%       large offsets.
%
%   (d) Mean f_ad per flight, which is the quantity Wood (2005) reports in his
%       Table 3. He notes that estimates from a single profile can be biased by
%       variations in cloud base (p. 3019), so a flight mean is the fairer
%       comparison with his values.
%
% INPUTS:
%   ensemble_profiles - the ensemble in any layout ensembleCell2StructOfCells
%                       accepts. Profiles need .adiabaticity (extrapolated base),
%                       .adiabaticity_measured_base and .adiabaticity_lowest_sample,
%                       written by find_verticalProfiles_VOCALS_REx_ver2.
%   campaign_label    - (optional) string for the title
%   rng_seed          - (optional, default 42) seed of the bootstrap. A private
%                       random stream is used, so the session's global random
%                       state is left untouched.
%
% OUTPUT:
%   out - structure with fields
%       .stats        table with one row per definition: median and percentiles
%                     of f_ad, the number above 1, the climb and descent medians,
%                     their offset with its 95% CI, and flights with mean f_ad > 1
%       .f_ad         one 1 x n_profiles vector per definition [-]
%       .dz_measured_m, .dz_extrapolate_m   base below the lowest sample [m]
%       .is_fallback                        extrapolation used its fallback
%       .is_ascending, .flight_id           per-profile labels
%       .rng_seed, .n_bootstrap             recorded for reproducibility
%
% SEE ALSO:
%   compute_adiabaticity_Wood2005, find_insitu_cloud_base, plot_ensemble_adiabaticity
%
% By Andrew John Buggee
%%

function out = plot_adiabaticity_base_comparison(ensemble_profiles, campaign_label, rng_seed)

if nargin < 2 || isempty(campaign_label)
    campaign_label = '';
end

if nargin < 3 || isempty(rng_seed)
    rng_seed = 42;
end

n_bootstrap = 10000;

%% Pull the per-profile numbers out of the ensemble

ep = ensembleCell2StructOfCells(ensemble_profiles);

required = {'adiabaticity', 'adiabaticity_measured_base', 'adiabaticity_lowest_sample'};
for kk = 1:numel(required)
    if ~isfield(ep, required{kk})
        error([newline, 'Profiles carry no .', required{kk}, ' field. Re-run ', ...
            'find_verticalProfiles_VOCALS_REx_ver2 to create it.', newline])
    end
end

n_prof = numel(ep.adiabaticity);

% The four versions compared. The last changes the numerator, not the base.
series_names = {'Lowest sample', 'Extrapolated base (default)', 'Measured base', 'Extrapolated base, droplet LWP'};
n_series = numel(series_names);

f_ad = nan(n_series, n_prof);          % [n_series x n_prof]
dz_measured_m = nan(1, n_prof);
dz_extrapolate_m = nan(1, n_prof);
is_fallback = false(1, n_prof);
is_ascending = false(1, n_prof);
flight_id = strings(1, n_prof);

for nn = 1:n_prof

    f_ad(1, nn) = ep.adiabaticity_lowest_sample{nn}.f_ad_lwp;
    f_ad(2, nn) = ep.adiabaticity{nn}.f_ad_lwp;
    f_ad(3, nn) = ep.adiabaticity_measured_base{nn}.f_ad_lwp;
    f_ad(4, nn) = ep.adiabaticity{nn}.f_ad_lwp_droplet;

    dz_measured_m(nn)    = ep.adiabaticity_measured_base{nn}.base_search.dz_below_original_m;
    dz_extrapolate_m(nn) = ep.adiabaticity{nn}.dz_base_m;
    is_fallback(nn)      = strcmp(ep.adiabaticity{nn}.base_source, 'fallback');
    is_ascending(nn)     = ep.adiabaticity_measured_base{nn}.base_search.is_ascending;
    flight_id(nn)        = string(ep.dateOfFlight{nn});

end


%% Statistics, including a bootstrap of the climb minus descent offset

% A private random stream, so results repeat and the global state is untouched
rand_stream = RandStream('mt19937ar', 'Seed', rng_seed);

idx_ascent  = find(is_ascending);
idx_descent = find(~is_ascending);

median_ascent  = nan(n_series, 1);
median_descent = nan(n_series, 1);
ci_ascent      = nan(n_series, 2);
ci_descent     = nan(n_series, 2);
offset_ci      = nan(n_series, 2);

% The same resamples are used for every definition, so they are compared fairly
resample_ascent  = randi(rand_stream, numel(idx_ascent),  numel(idx_ascent),  n_bootstrap);   % [n_ascent x n_bootstrap]
resample_descent = randi(rand_stream, numel(idx_descent), numel(idx_descent), n_bootstrap);   % [n_descent x n_bootstrap]

for ss = 1:n_series

    f_up   = f_ad(ss, idx_ascent);
    f_down = f_ad(ss, idx_descent);

    median_ascent(ss)  = median(f_up,   'omitnan');
    median_descent(ss) = median(f_down, 'omitnan');

    boot_up   = median(reshape(f_up(resample_ascent),   size(resample_ascent)),  1, 'omitnan');
    boot_down = median(reshape(f_down(resample_descent), size(resample_descent)), 1, 'omitnan');

    ci_ascent(ss, :)  = prctile(boot_up,   [2.5, 97.5]);
    ci_descent(ss, :) = prctile(boot_down, [2.5, 97.5]);
    offset_ci(ss, :)  = prctile(boot_up - boot_down, [2.5, 97.5]);

end

% Mean f_ad per flight, which is what Wood (2005) reports in his Table 3
[flights, ~, flight_index] = unique(flight_id, 'stable');
n_flights = numel(flights);
f_flight = nan(n_series, n_flights);
for ff = 1:n_flights
    f_flight(:, ff) = mean(f_ad(:, flight_index == ff), 2, 'omitnan');
end

definition      = string(series_names(:));
median_f_ad     = median(f_ad, 2, 'omitnan');
p10_f_ad        = prctile(f_ad, 10, 2);
p90_f_ad        = prctile(f_ad, 90, 2);
n_above_1       = sum(f_ad > 1, 2);
offset          = median_ascent - median_descent;
offset_ci_low   = offset_ci(:, 1);
offset_ci_high  = offset_ci(:, 2);
flights_above_1 = sum(f_flight > 1, 2);

stats = table(definition, median_f_ad, p10_f_ad, p90_f_ad, n_above_1, median_ascent, ...
    median_descent, offset, offset_ci_low, offset_ci_high, flights_above_1);

out.stats = stats;
out.f_ad = f_ad;
out.series_names = series_names;
out.dz_measured_m = dz_measured_m;
out.dz_extrapolate_m = dz_extrapolate_m;
out.is_fallback = is_fallback;
out.is_ascending = is_ascending;
out.flight_id = flight_id;
out.f_flight = f_flight;
out.flights = flights;
out.rng_seed = rng_seed;
out.n_bootstrap = n_bootstrap;

fprintf('\nCloud-base comparison: %d profiles (%d climbs, %d descents), %d flights. ', ...
    n_prof, numel(idx_ascent), numel(idx_descent), n_flights);
fprintf('Bootstrap: B = %d, seed = %d\n', n_bootstrap, rng_seed);
disp(stats)


%% Figure

% Categorical colors in a fixed order, one per cloud-base definition. The
% droplet-LWP series keeps the color of the extrapolated base because it uses the
% same base; line style and marker set it apart.
clr = [0.165, 0.471, 0.839;     % blue    #2a78d6  lowest sample
       0.922, 0.408, 0.204;     % orange  #eb6834  extrapolated base
       0.106, 0.686, 0.478;     % aqua    #1baf7a  measured base
       0.922, 0.408, 0.204];    % orange  #eb6834  extrapolated base, droplet LWP
line_style = {'-', '-', '-', '--'};
marker     = {'o', 'o', 'o', 'd'};
clr_reference = [0.45, 0.45, 0.45];
clr_points    = [0.25, 0.25, 0.25];

fig = new_factory_figure('Position', [50, 50, 1250, 850]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


% ---- (a) Cumulative distribution of per-profile f_ad --------------------
ax = nexttile(tl);
hold(ax, 'on')

x_max = max(2, min(3, prctile(f_ad(:), 99)));
h = gobjects(n_series, 1);

for ss = 1:n_series
    f_sorted = sort(f_ad(ss, isfinite(f_ad(ss, :))));
    cdf = (1:numel(f_sorted)) ./ numel(f_sorted);
    h(ss) = stairs(ax, f_sorted, cdf, 'Color', clr(ss, :), 'LineStyle', line_style{ss}, 'LineWidth', 2);
end

xline(ax, 1, 'Color', clr_reference, 'LineWidth', 1)

legend_str = cell(n_series, 1);
for ss = 1:n_series
    legend_str{ss} = sprintf('%s: median %.2f, %d above 1', series_names{ss}, median_f_ad(ss), n_above_1(ss));
end
legend(ax, h, legend_str, 'Location', 'southeast', 'Box', 'off', 'FontSize', 9)

xlim(ax, [0, x_max]); ylim(ax, [0, 1])
xlabel(ax, 'f_{ad} = LWP_{obs} / LWP_{ad}')
ylabel(ax, 'Fraction of profiles at or below f_{ad}')
title(ax, '(a) Per-profile adiabaticity', 'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1;


% ---- (b) Where each definition places the base ---------------------------
ax = nexttile(tl);
hold(ax, 'on')

axis_low  = min([dz_measured_m, dz_extrapolate_m, 0]);
axis_high = max([dz_measured_m, dz_extrapolate_m, 50]);
axis_pad  = 0.05 * (axis_high - axis_low);      % keep points off the frame
axis_low  = axis_low - axis_pad;
axis_high = axis_high + axis_pad;

plot(ax, [axis_low, axis_high], [axis_low, axis_high], '-', 'Color', clr_reference, 'LineWidth', 1, ...
    'HandleVisibility', 'off')

scatter(ax, dz_extrapolate_m(is_ascending), dz_measured_m(is_ascending), 36, clr_points, 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Climb')
scatter(ax, dz_extrapolate_m(~is_ascending), dz_measured_m(~is_ascending), 36, clr_points, ...
    'LineWidth', 1, 'DisplayName', 'Descent')
scatter(ax, dz_extrapolate_m(is_fallback), dz_measured_m(is_fallback), 90, clr(2, :), 'x', ...
    'LineWidth', 1.5, 'DisplayName', sprintf('Fallback used (%d)', sum(is_fallback)))

legend(ax, 'Location', 'northwest', 'Box', 'off')
xlim(ax, [axis_low, axis_high]); ylim(ax, [axis_low, axis_high])
xlabel(ax, 'Extrapolated base: distance below lowest sample (m)')
ylabel(ax, 'Measured base: distance below lowest sample (m)')
title(ax, sprintf('(b) Medians: measured %.0f m, extrapolated %.0f m', ...
    median(dz_measured_m, 'omitnan'), median(dz_extrapolate_m, 'omitnan')), 'FontWeight', 'normal')
text(ax, 0.97, 0.06, 'Negative: profile trimmed, it ran below the base', 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', [0.35, 0.35, 0.35])
grid(ax, 'on'); ax.GridAlpha = 0.1;


% ---- (c) Climbs versus descents: the test of each definition --------------
ax = nexttile(tl);
hold(ax, 'on')

row_gap = 0.17;     % vertical separation of the climb and descent markers of a row
y_row = n_series:-1:1;

for ss = 1:n_series

    % Climb: filled marker, above the row center
    plot(ax, ci_ascent(ss, :), [y_row(ss), y_row(ss)] + row_gap, '-', 'Color', clr(ss, :), 'LineWidth', 2, ...
        'LineStyle', line_style{ss})
    plot(ax, median_ascent(ss), y_row(ss) + row_gap, marker{ss}, 'MarkerSize', 8, ...
        'MarkerFaceColor', clr(ss, :), 'MarkerEdgeColor', clr(ss, :))

    % Descent: open marker, below the row center
    plot(ax, ci_descent(ss, :), [y_row(ss), y_row(ss)] - row_gap, '-', 'Color', clr(ss, :), 'LineWidth', 2, ...
        'LineStyle', line_style{ss})
    plot(ax, median_descent(ss), y_row(ss) - row_gap, marker{ss}, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', clr(ss, :), 'LineWidth', 1.5)

end

xline(ax, 1, 'Color', clr_reference, 'LineWidth', 1)

% The offset is what matters here, so it is written beside each row
row_labels = cell(n_series, 1);
for ss = 1:n_series
    row_labels{n_series - ss + 1} = sprintf('%s\\newline\\Delta = %+.2f [%+.2f, %+.2f]', ...
        series_names{ss}, offset(ss), offset_ci_low(ss), offset_ci_high(ss));
end

yticks(ax, 1:n_series); yticklabels(ax, row_labels);
ylim(ax, [0.4, n_series + 0.6])
xlabel(ax, 'Median f_{ad} with bootstrap 95% CI')
title(ax, '(c) Climbs (filled) vs descents (open); \Delta = climb - descent', 'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1; ax.YGrid = 'off';


% ---- (d) Flight means, the quantity Wood (2005) reports -------------------
ax = nexttile(tl);
hold(ax, 'on')

dodge = linspace(-0.24, 0.24, n_series);
h = gobjects(n_series, 1);

for ss = 1:n_series
    is_open = ss == 4;
    face = clr(ss, :);
    if is_open
        face = 'w';
    end
    h(ss) = plot(ax, (1:n_flights) + dodge(ss), f_flight(ss, :), marker{ss}, 'MarkerSize', 7, ...
        'MarkerFaceColor', face, 'MarkerEdgeColor', clr(ss, :), 'LineWidth', 1.2);
end

yline(ax, 1, 'Color', clr_reference, 'LineWidth', 1)

flight_labels = extractBefore(flights, 7);      % '15-Oct-2008' -> '15-Oct'
xticks(ax, 1:n_flights); xticklabels(ax, flight_labels); xtickangle(ax, 45)
xlim(ax, [0.4, n_flights + 0.6])
ylim(ax, [0, max(1.5, 1.05 * max(f_flight(:)))])
ylabel(ax, 'Flight mean f_{ad}')
legend(ax, h, series_names, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9, 'NumColumns', 2)
title(ax, sprintf('(d) Flight means; above 1: %d, %d, %d, %d of %d flights', flights_above_1, n_flights), ...
    'FontWeight', 'normal')
grid(ax, 'on'); ax.GridAlpha = 0.1; ax.XGrid = 'off';


if isempty(campaign_label)
    title(tl, sprintf('Cloud-base definitions for f_{ad}: %d profiles', n_prof))
else
    title(tl, sprintf('Cloud-base definitions for f_{ad}: %s, %d profiles', campaign_label, n_prof))
end

out.figure = fig;

end
