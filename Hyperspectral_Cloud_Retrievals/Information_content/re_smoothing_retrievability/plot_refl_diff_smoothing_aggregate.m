function plot_refl_diff_smoothing_aggregate(results_dir, varargin)
%% Aggregate retrievability figure across all in-situ profiles
%
% Loads every per-profile result produced by hysics_refl_smoothingTest_from_insitu
% (file pattern: refl_smoothingTest_*.mat) in results_dir, then builds the
% summary figure: one panel per smoothing window showing, as a function of
% wavelength, the distribution across profiles of the reflectance difference
% |R_raw - R_smoothed| together with the HySICS measurement-uncertainty
% envelope (0.3% of the reflectance).
%
% Interpretation: wherever the |dR| distribution sits ABOVE the uncertainty
% envelope, smoothing that much vertical droplet-size structure changes the
% HySICS reflectance by more than the instrument noise -- i.e. the vertical
% structure is detectable (a sufficient condition for retrievability). Where
% it sits below, those droplet-profile differences are NOT retrievable.
%
% INPUTS:
%   results_dir - folder containing the refl_smoothingTest_*.mat files.
%                 If omitted, defaults to a local results folder.
%
% Name-value options:
%   'SaveFig'   - logical, save the figure as a .png (default false)
%   'YScale'    - 'linear' (default) or 'log'
%   'Band'      - percentile band half-spans to shade, default [25 75] and
%                 [5 95] are both drawn. Pass [] to draw only the median.
%
% Example:
%   plot_refl_diff_smoothing_aggregate('/path/to/results/')
%
% By Andrew John Buggee

%% ----------------------------- parse inputs -------------------------------

p = inputParser;
addRequired(p,  'results_dir', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveFig', false, @islogical);
addParameter(p, 'YScale', 'linear', @(s) any(strcmpi(s, {'linear','log'})));
addParameter(p, 'ShowBands', true, @islogical);
parse(p, results_dir, varargin{:});
results_dir = char(p.Results.results_dir);
if results_dir(end) ~= filesep; results_dir = [results_dir, filesep]; end

%% ----------------------------- load results -------------------------------

files = dir([results_dir, 'refl_smoothingTest_*.mat']);
if isempty(files)
    error([newline, 'No refl_smoothingTest_*.mat files found in: ', results_dir, newline])
end
fprintf('Found %d result files in %s\n', numel(files), results_dir);

% load the first to establish the common grid
first = load([results_dir, files(1).name], 'wl_nm', 'smoothing_windows_m');
wl_nm    = first.wl_nm(:);
windows_m = first.smoothing_windows_m(:)';
num_wl  = numel(wl_nm);
num_win = numel(windows_m);

% accumulate |dR| and the uncertainty envelope across profiles
absdR_all   = nan(num_wl, num_win, numel(files));   % |R_raw - R_smoothed_w|
uncert_raw  = nan(num_wl, numel(files));            % 0.003 * R_raw
campaign_of = strings(numel(files), 1);
n_used = 0;

for ff = 1:numel(files)
    S = load([results_dir, files(ff).name], ...
        'wl_nm', 'dR', 'Refl_uncert_hysics', 'smoothing_windows_m', 'campaign_name');

    % skip files that don't match the common grid / windows
    if numel(S.wl_nm) ~= num_wl || numel(S.smoothing_windows_m) ~= num_win || ...
            any(abs(S.smoothing_windows_m(:)' - windows_m) > 1e-9)
        warning('Skipping %s (grid/window mismatch).', files(ff).name);
        continue
    end

    n_used = n_used + 1;
    absdR_all(:, :, n_used) = abs(S.dR);
    uncert_raw(:, n_used)   = S.Refl_uncert_hysics(:, 1);   % envelope of the raw spectrum
    campaign_of(n_used)     = string(S.campaign_name);
end

absdR_all  = absdR_all(:, :, 1:n_used);
uncert_raw = uncert_raw(:, 1:n_used);
fprintf('Aggregated %d profiles (%d VOCALS, %d ORACLES).\n', n_used, ...
    sum(campaign_of(1:n_used)=="vocalsRex"), sum(campaign_of(1:n_used)=="oracles"));

%% --------------------------- summary statistics ---------------------------

med_absdR = median(absdR_all, 3, 'omitnan');                 % num_wl x num_win
p25 = prctile(absdR_all, 25, 3);
p75 = prctile(absdR_all, 75, 3);
p05 = prctile(absdR_all,  5, 3);
p95 = prctile(absdR_all, 95, 3);

med_uncert = median(uncert_raw, 2, 'omitnan');               % num_wl x 1

% fraction of profiles whose |dR| exceeds their own uncertainty, per band
uncert_3d = reshape(uncert_raw, num_wl, 1, n_used);
frac_retrievable = mean(absdR_all > uncert_3d, 3, 'omitnan'); % num_wl x num_win

%% --------------------------------- plot -----------------------------------

figure('Name', 'r_e smoothing retrievability', ...
    'Position', [60, 60, 1500, 420]);
tl = tiledlayout(1, num_win, 'TileSpacing', 'compact', 'Padding', 'compact');

band_color   = [0.30 0.45 0.80];
median_color = [0.10 0.20 0.55];

for ww = 1:num_win
    nexttile; hold on

    if p.Results.ShowBands
        % 5-95th percentile band (light)
        fill([wl_nm; flipud(wl_nm)], [p05(:,ww); flipud(p95(:,ww))], band_color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        % 25-75th percentile band (darker)
        fill([wl_nm; flipud(wl_nm)], [p25(:,ww); flipud(p75(:,ww))], band_color, ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none', 'DisplayName', '25-75th pctile');
    end

    % median |dR|
    plot(wl_nm, med_absdR(:,ww), '-', 'Color', median_color, 'LineWidth', 2, ...
        'DisplayName', 'median |\DeltaR|');

    % HySICS uncertainty envelope (median across profiles)
    plot(wl_nm, med_uncert, '--', 'Color', [0.80 0.10 0.10], 'LineWidth', 2, ...
        'DisplayName', 'HySICS 0.3% uncert.');

    grid on; box on
    set(gca, 'YScale', p.Results.YScale)
    xlim([min(wl_nm), max(wl_nm)])
    xlabel('Wavelength  [nm]')
    if ww == 1; ylabel('|R_{raw} - R_{smoothed}|'); end
    title(sprintf('window = %d m', windows_m(ww)))
    if ww == num_win; legend('Location', 'northwest'); end
end

title(tl, sprintf(['HySICS reflectance sensitivity to vertical r_e smoothing ', ...
    '(%d in-situ profiles, nadir)'], n_used), 'FontWeight', 'bold')

if p.Results.SaveFig
    saveas(gcf, [results_dir, 'aggregate_refl_diff_vs_smoothing.png']);
    fprintf('Saved figure to %saggregate_refl_diff_vs_smoothing.png\n', results_dir);
end

%% ------------------- print a compact retrievability summary ----------------

fprintf('\n==== Retrievability summary (median across %d profiles) ====\n', n_used);
fprintf('  window [m] | %% of HySICS bands where median|dR| > 0.3%% uncert\n');
for ww = 1:num_win
    pct_bands = 100 * mean(med_absdR(:,ww) > med_uncert);
    fprintf('  %8d   |                 %5.1f%%\n', windows_m(ww), pct_bands);
end

% Optional companion figure: fraction of profiles retrievable per band
figure('Name', 'Fraction of profiles retrievable', 'Position', [80, 520, 1500, 360]);
tl2 = tiledlayout(1, num_win, 'TileSpacing', 'compact', 'Padding', 'compact');
for ww = 1:num_win
    nexttile
    plot(wl_nm, 100*frac_retrievable(:,ww), '-', 'Color', median_color, 'LineWidth', 1.5);
    yline(50, ':k');
    grid on; box on; ylim([0 100]); xlim([min(wl_nm), max(wl_nm)])
    xlabel('Wavelength  [nm]')
    if ww == 1; ylabel('% of profiles with |\DeltaR| > uncert'); end
    title(sprintf('window = %d m', windows_m(ww)))
end
title(tl2, 'Fraction of profiles whose smoothing is detectable, per band', ...
    'FontWeight', 'bold')

if p.Results.SaveFig
    saveas(gcf, [results_dir, 'aggregate_fraction_retrievable.png']);
end

end
