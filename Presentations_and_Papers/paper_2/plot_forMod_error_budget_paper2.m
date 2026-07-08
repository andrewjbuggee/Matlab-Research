%% Analyze and plot the forward-model-parameter error budget for the EMIT-Aqua retrievals
%
% This is the *analysis / visualization* companion to
% compute_forMod_error_budget_paper2.m. It takes the per-pixel forward-model
% error decomposition produced there and views it several ways so the paper-2
% error budget can be interpreted:
%
%   FIG 1  Error-budget bar chart
%          Median fractional 1-sigma uncertainty contributed by each physical
%          forward-model source (r_e profile, cloud top height, effective
%          variance), per retrieved variable, with 16th/84th-percentile whiskers
%          and the forward-model-total and posterior-total overlaid.
%          -> "How big is each source, in relative terms?"
%
%   FIG 2  Share-of-variance (the statistically correct "fraction of total")
%          (a) 100%-stacked: each source as a share of the FORWARD-MODEL variance.
%          (b) forward-model total as a share of the TOTAL posterior variance.
%          Shares are computed in VARIANCE (sigma^2) because uncertainties add in
%          quadrature: sigma_total^2 = sigma_meas^2 + sigma_fm^2 + ... , so the
%          additive share of a source is sigma_source^2 / sigma_total^2. Doing
%          this with sigma directly would overcount.
%          -> "What fraction of the retrieval uncertainty is forward model, and
%              which source dominates it?"
%
%   FIG 3  Uncertainty vs retrieved optical thickness tau_c (quantile-binned)
%          Median (+ shaded 16-84 band) of each source's contribution as a
%          function of retrieved tau_c, one subplot per retrieved variable.
%          Shown in both fractional and absolute units.
%          -> "Does the forward-model error grow or shrink with tau_c?"
%
%   FIG 4  Forward-model share of total variance vs retrieved tau_c
%          One line per retrieved variable.
%          -> "At what tau_c does the forward model start to matter?"
%
% WHY tau_c is binned by QUANTILE, not fixed width: retrieved tau_c is skewed,
% so equal-count (quantile) bins keep the per-bin median equally well determined
% across the range instead of leaving sparse high-tau bins.
%
% KEY STRUCTURAL FACT (from calc_retrieval_..._perPixel.m, line ~1101):
%   total_meas_fm_cov = measurement_cov_log + K_b * S_b * K_b'
% so the forward-model error is folded INTO the effective measurement covariance
% and posterior_cov_log already contains it. The forward-model term
%   S_f = G * K_b * S_b * K_b' * G'
% is therefore a genuine sub-component of the total posterior covariance, which
% is what makes the FIG-2 / FIG-4 "fraction of total" ratios meaningful.
%
% Data source: prefers the newest EMIT_Aqua_forMod_error_budget_*.mat saved by
% compute_forMod_error_budget_paper2.m. If none is found it runs the identical
% per-pixel extraction itself so this script is self-contained.
%
% By Andrew John Buggee
%%

clear variables

% --------------------------- configuration ----------------------------
which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374') == true
    % ------ Mac Desktop ------
    retrieval_directory = '/Users/anbu8374/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_13/';
    save_directory      = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/';

elseif strcmp(which_computer, 'andrewbuggee') == true
    % ------ Macbook ------
    retrieval_directory = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_13/';
    save_directory      = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/';

else
    error([newline, 'Define retrieval_directory / save_directory for this computer.', newline])
end

n_tau_bins   = 6;        % number of quantile (equal-count) bins along retrieved tau_c
min_bin_count = 5;       % drop bins with fewer than this many pixels
save_figures = false;    % set true to write PNG/FIG into save_directory/saved_figures
% ----------------------------------------------------------------------


% ---- variable and source bookkeeping (must match the compute script) ----
% Row order of the retrieval / posterior_cov_log: ln(r_top), ln(r_bot), ln(tau_c), ln(acpw)
var_fields = {'r_top', 'r_bot', 'tau_c', 'acpw'};
var_labels = {'r_{top}', 'r_{bot}', '\tau_c', 'acpw'};
var_units  = {'\mum', '\mum', '', 'cm'};      % absolute-unit labels for retrieved values
n_vars     = numel(var_fields);
i_tau      = find(strcmp(var_fields, 'tau_c'));   % index of tau_c (the x-axis variable)

source_fields = {'re_profile', 'cloud_top_height', 'effective_variance'};
source_labels = {'r_e profile (adiabatic)', 'cloud top height', 'effective variance'};
n_src         = numel(source_fields);

% consistent colors: 3 sources, then FM-total (black), posterior-total (red)
src_colors = [0.20 0.45 0.70;    % blue   - r_e profile
              0.90 0.60 0.10;    % orange - cloud top height
              0.35 0.70 0.35];   % green  - effective variance
col_fm_total  = [0 0 0];
col_tot_post  = [0.85 0.15 0.15];


%% ---- get the per-pixel arrays (load saved, else extract) ----
[frac_total, frac_fm_total, frac_by_source, ...
 abs_total,  abs_fm_total,  abs_by_source, xhat_all, n_used, src_from_data] = ...
    get_perpixel_arrays(save_directory, retrieval_directory, source_fields, n_vars, n_src);

% guard: the loaded/extracted source order must match what this script assumes
if ~isempty(src_from_data) && ~isequal(src_from_data(:)', source_fields)
    warning(['Source-field order in the data (%s) differs from this script''s ', ...
        'assumption (%s). Check source_fields.'], strjoin(src_from_data, ', '), strjoin(source_fields, ', '));
end

% ---- keep only fully finite pixels ----
valid = all(isfinite(xhat_all), 1) & all(isfinite(frac_total), 1) & all(isfinite(frac_fm_total), 1);
frac_total     = frac_total(:, valid);
frac_fm_total  = frac_fm_total(:, valid);
frac_by_source = frac_by_source(:, :, valid);
abs_total      = abs_total(:, valid);
abs_fm_total   = abs_fm_total(:, valid);
abs_by_source  = abs_by_source(:, :, valid);
xhat_all       = xhat_all(:, valid);
N              = size(xhat_all, 2);

fprintf('\nAnalyzing %d valid pixels (of %d used in extraction).\n', N, n_used);
if N < min_bin_count
    warning('Only %d valid pixels; tau-binned figures (3,4) may be unreliable.', N);
end


%% ---- derived per-pixel quantities (all in VARIANCE space) ----
% variance (sigma^2) versions
var_total     = frac_total.^2;                 % n_vars x N
var_fm_total  = frac_fm_total.^2;              % n_vars x N
var_by_source = frac_by_source.^2;             % n_vars x n_src x N

% share of each source within the forward-model budget (sums to 1 over sources)
share_src_of_fm = var_by_source ./ max(sum(var_by_source, 2), eps);   % n_vars x n_src x N

% forward-model total as a share of the TOTAL posterior variance
share_fm_of_total = var_fm_total ./ max(var_total, eps);              % n_vars x N

tau_c_retrieved = xhat_all(i_tau, :);          % 1 x N, retrieved optical thickness


%% =====================================================================
%  FIGURE 1 — error-budget bar chart (fractional 1-sigma, %)
%  =====================================================================
med_src_frac = 100 * median(frac_by_source, 3, 'omitnan');           % n_vars x n_src (%)
p_src_frac   = 100 * prctile(frac_by_source, [16 84], 3);            % n_vars x n_src x 2
med_fm_frac  = 100 * median(frac_fm_total, 2, 'omitnan');            % n_vars x 1
med_tot_frac = 100 * median(frac_total,    2, 'omitnan');            % n_vars x 1

figure('Position', [80 80 900 520]);
b = bar(med_src_frac, 'grouped');
for s = 1:n_src
    b(s).FaceColor = src_colors(s, :);
    b(s).DisplayName = source_labels{s};
end
hold on

% 16-84 whiskers on each grouped bar
for s = 1:n_src
    xc = b(s).XEndPoints;                                   % x of each bar in this group
    lo = 100 * p_src_frac(:, s, 1).';
    hi = 100 * p_src_frac(:, s, 2).';
    errorbar(xc, med_src_frac(:, s).', med_src_frac(:, s).' - lo, hi - med_src_frac(:, s).', ...
        'LineStyle', 'none', 'Color', 0.3*[1 1 1], 'CapSize', 4, 'HandleVisibility', 'off');
end

h_fm  = plot(1:n_vars, med_fm_frac,  's', 'MarkerFaceColor', col_fm_total, 'MarkerEdgeColor', col_fm_total, ...
    'MarkerSize', 9, 'LineStyle', 'none', 'DisplayName', 'forward-model (all sources)');
h_tot = plot(1:n_vars, med_tot_frac, 'd', 'MarkerFaceColor', col_tot_post, 'MarkerEdgeColor', col_tot_post, ...
    'MarkerSize', 9, 'LineStyle', 'none', 'DisplayName', 'total posterior');

set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Fractional 1\sigma uncertainty (%)', 'Interpreter', 'tex');
title(sprintf('EMIT-Aqua forward-model error budget (median of %d pixels; whiskers 16-84th pct)', N));
legend([b, h_fm, h_tot], 'Location', 'best', 'Interpreter', 'tex');
grid on; grid minor
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig1_bar');


%% =====================================================================
%  FIGURE 2 — share of variance (the correct "fraction of total")
%  =====================================================================
figure('Position', [110 110 1000 460]);

% --- (a) each source as a share of the forward-model variance (sums to 100%) ---
subplot(1, 2, 1);
med_share_src_of_fm = 100 * median(share_src_of_fm, 3, 'omitnan');    % n_vars x n_src (%)
bs = bar(med_share_src_of_fm, 'stacked');
for s = 1:n_src
    bs(s).FaceColor = src_colors(s, :);
    bs(s).DisplayName = source_labels{s};
end
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Share of forward-model variance (%)', 'Interpreter', 'tex');
ylim([0 100]);
title('(a) Which source dominates the forward-model budget');
legend(bs, 'Location', 'southoutside', 'Interpreter', 'tex', 'NumColumns', 1);
grid on

% --- (b) forward-model total as a share of the TOTAL posterior variance ---
subplot(1, 2, 2);
med_share_fm_of_total = 100 * median(share_fm_of_total, 2, 'omitnan');   % n_vars x 1
p_share_fm_of_total   = 100 * prctile(share_fm_of_total, [16 84], 2);    % n_vars x 2
bb = bar(1:n_vars, med_share_fm_of_total, 0.6, 'FaceColor', [0.5 0.5 0.5]);
hold on
errorbar(1:n_vars, med_share_fm_of_total, ...
    med_share_fm_of_total - p_share_fm_of_total(:, 1), p_share_fm_of_total(:, 2) - med_share_fm_of_total, ...
    'LineStyle', 'none', 'Color', 'k', 'CapSize', 6);
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Forward-model share of total variance (%)', 'Interpreter', 'tex');
title('(b) How much of the total retrieval error is forward model');
grid on
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig2_shares');


%% =====================================================================
%  FIGURE 3 — uncertainty vs retrieved tau_c (quantile-binned)
%  Two rows: top = fractional (%), bottom = absolute (retrieved units)
%  =====================================================================
[bin_ctr, bin_lo, bin_hi, bin_idx] = quantile_bins(tau_c_retrieved, n_tau_bins, min_bin_count);
n_bins_kept = numel(bin_ctr);

figure('Position', [140 60 1200 720]);
for i = 1:n_vars

    % ---------- top row: fractional (%) ----------
    subplot(2, n_vars, i);
    hold on
    lg = gobjects(1, n_src + 2);
    for s = 1:n_src
        y = squeeze(frac_by_source(i, s, :)).' * 100;       % 1 x N (%)
        [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
        shade_band(bin_ctr, plo, phi, src_colors(s, :));
        lg(s) = plot(bin_ctr, m, '-o', 'Color', src_colors(s, :), 'MarkerFaceColor', src_colors(s, :), ...
            'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', source_labels{s});
    end
    % forward-model total and posterior total medians
    [m_fm, ~, ~]   = bin_stats(frac_fm_total(i, :) * 100, bin_idx, n_bins_kept);
    [m_tot, ~, ~]  = bin_stats(frac_total(i, :)    * 100, bin_idx, n_bins_kept);
    lg(n_src+1) = plot(bin_ctr, m_fm,  's--', 'Color', col_fm_total, 'MarkerFaceColor', col_fm_total, ...
        'MarkerSize', 4, 'DisplayName', 'forward-model (all)');
    lg(n_src+2) = plot(bin_ctr, m_tot, 'd--', 'Color', col_tot_post, 'MarkerFaceColor', col_tot_post, ...
        'MarkerSize', 4, 'DisplayName', 'total posterior');
    xlabel('retrieved \tau_c', 'Interpreter', 'tex');
    ylabel('fractional 1\sigma (%)', 'Interpreter', 'tex');
    title(var_labels{i}, 'Interpreter', 'tex');
    grid on
    if i == 1, legend(lg, 'Location', 'best', 'Interpreter', 'tex', 'FontSize', 7); end

    % ---------- bottom row: absolute (retrieved units) ----------
    subplot(2, n_vars, n_vars + i);
    hold on
    for s = 1:n_src
        y = squeeze(abs_by_source(i, s, :)).';              % 1 x N (absolute)
        [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
        shade_band(bin_ctr, plo, phi, src_colors(s, :));
        plot(bin_ctr, m, '-o', 'Color', src_colors(s, :), 'MarkerFaceColor', src_colors(s, :), ...
            'MarkerSize', 4, 'LineWidth', 1.5);
    end
    [m_fm_a, ~, ~]  = bin_stats(abs_fm_total(i, :), bin_idx, n_bins_kept);
    [m_tot_a, ~, ~] = bin_stats(abs_total(i, :),    bin_idx, n_bins_kept);
    plot(bin_ctr, m_fm_a,  's--', 'Color', col_fm_total, 'MarkerFaceColor', col_fm_total, 'MarkerSize', 4);
    plot(bin_ctr, m_tot_a, 'd--', 'Color', col_tot_post, 'MarkerFaceColor', col_tot_post, 'MarkerSize', 4);
    xlabel('retrieved \tau_c', 'Interpreter', 'tex');
    if isempty(var_units{i})
        ylabel('absolute 1\sigma', 'Interpreter', 'tex');
    else
        ylabel(sprintf('absolute 1\\sigma (%s)', var_units{i}), 'Interpreter', 'tex');
    end
    grid on
end
sgtitle(sprintf('Uncertainty vs retrieved \\tau_c  (%d quantile bins; shaded = 16-84th pct)', n_bins_kept), ...
    'Interpreter', 'tex');
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig3_vsTau');


%% =====================================================================
%  FIGURE 4 — forward-model share of total variance vs retrieved tau_c
%  =====================================================================
figure('Position', [170 90 780 520]);
hold on
line_cols = lines(n_vars);
for i = 1:n_vars
    y = share_fm_of_total(i, :) * 100;                       % 1 x N (%)
    [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
    shade_band(bin_ctr, plo, phi, line_cols(i, :));
    plot(bin_ctr, m, '-o', 'Color', line_cols(i, :), 'MarkerFaceColor', line_cols(i, :), ...
        'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName', var_labels{i});
end
xlabel('retrieved \tau_c', 'Interpreter', 'tex');
ylabel('forward-model share of total variance (%)', 'Interpreter', 'tex');
title(sprintf('When does the forward model matter? (%d quantile bins; shaded 16-84th pct)', n_bins_kept));
legend('Location', 'best', 'Interpreter', 'tex');
grid on; grid minor
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig4_shareVsTau');


%% ---- concise text summary ----
fprintf('\n==================== forward-model error budget (median of %d pixels) ====================\n', N);
for i = 1:n_vars
    fprintf('\n%-6s (retrieved %s median = %.3g)\n', var_fields{i}, var_labels{i}, median(xhat_all(i, :), 'omitnan'));
    fprintf('   forward-model 1sigma : %6.2f%% of value   |   %6.1f%% of total variance\n', ...
        med_fm_frac(i), 100*median(share_fm_of_total(i, :), 'omitnan'));
    for s = 1:n_src
        fprintf('     - %-22s : %6.2f%% of value   |   %6.1f%% of FM variance\n', ...
            source_labels{s}, med_src_frac(i, s), 100*median(share_src_of_fm(i, s, :), 'omitnan'));
    end
end
fprintf('\n(Shares are variance-based; sources sum to ~100%% of the forward-model variance.)\n\n');


%% ============================ local functions ============================

function [frac_total, frac_fm_total, frac_by_source, abs_total, abs_fm_total, ...
          abs_by_source, xhat_all, n_used, src_from_data] = ...
          get_perpixel_arrays(save_directory, retrieval_directory, source_fields, n_vars, n_src)
% Prefer the newest saved budget file; otherwise re-extract from the retrievals.

    src_from_data = {};

    saved = dir(fullfile(save_directory, 'EMIT_Aqua_forMod_error_budget_*.mat'));
    if ~isempty(saved)
        [~, newest] = max([saved.datenum]);
        f = fullfile(saved(newest).folder, saved(newest).name);
        L = load(f, 'frac_total', 'frac_fm_total', 'frac_by_source', ...
            'abs_total', 'abs_fm_total', 'abs_by_source', 'xhat_all', 'n_used', 'budget');
        needed = {'frac_total', 'frac_fm_total', 'frac_by_source', 'abs_total', ...
            'abs_fm_total', 'abs_by_source', 'xhat_all'};
        if all(isfield(L, needed))
            fprintf('Loaded per-pixel arrays from saved budget file:\n  %s\n', f);
            frac_total    = L.frac_total;     frac_fm_total = L.frac_fm_total;
            frac_by_source = L.frac_by_source; abs_total    = L.abs_total;
            abs_fm_total  = L.abs_fm_total;   abs_by_source = L.abs_by_source;
            xhat_all      = L.xhat_all;
            if isfield(L, 'n_used'), n_used = L.n_used; else, n_used = size(xhat_all, 2); end
            if isfield(L, 'budget') && isfield(L.budget, 'source_fields')
                src_from_data = L.budget.source_fields;
            end
            return
        end
        warning('Saved budget file is missing per-pixel arrays; re-extracting from retrievals.');
    end

    % ---- fallback: identical extraction to compute_forMod_error_budget_paper2.m ----
    fprintf('No usable saved budget file; extracting per-pixel arrays from:\n  %s\n', retrieval_directory);
    file_list = dir(fullfile(retrieval_directory, 'EMIT_dropRetrieval_*.mat'));
    if isempty(file_list)
        error([newline, 'No EMIT_dropRetrieval_*.mat files found in:', newline, retrieval_directory, newline]);
    end
    n_files = numel(file_list);

    frac_total     = nan(n_vars, n_files);
    frac_fm_total  = nan(n_vars, n_files);
    frac_by_source = nan(n_vars, n_src, n_files);
    abs_total      = nan(n_vars, n_files);
    abs_fm_total   = nan(n_vars, n_files);
    abs_by_source  = nan(n_vars, n_src, n_files);
    xhat_all       = nan(n_vars, n_files);
    n_used = 0;

    needed = {'retrieval', 'posterior_cov_log', 'forward_model_parameter_cov', ...
        'forMod_sensitivity_log', 'forMod_param_variance_log', 'forMod_param_names'};

    for k = 1:n_files
        S = load(fullfile(file_list(k).folder, file_list(k).name), 'GN_outputs');
        if ~isfield(S, 'GN_outputs') || isempty(S.GN_outputs), continue; end
        GN = S.GN_outputs;
        if ~all(isfield(GN, needed)), continue; end

        x_hat = GN.retrieval(:, end);
        if numel(x_hat) ~= n_vars || any(~isfinite(x_hat)), continue; end

        sigma_frac_total = sqrt(diag(GN.posterior_cov_log));
        frac_total(:, k) = sigma_frac_total;
        abs_total(:, k)  = x_hat .* sigma_frac_total;
        xhat_all(:, k)   = x_hat;

        sigma_frac_fm = sqrt(diag(GN.forward_model_parameter_cov));
        frac_fm_total(:, k) = sigma_frac_fm;
        abs_fm_total(:, k)  = x_hat .* sigma_frac_fm;

        sens  = GN.forMod_sensitivity_log;
        var_b = GN.forMod_param_variance_log(:).';
        names = GN.forMod_param_names;
        contrib_var = (sens.^2) .* var_b;

        for s = 1:n_src
            cols = strcmp(names, source_fields{s});
            var_s = sum(contrib_var(:, cols), 2);
            frac_by_source(:, s, k) = sqrt(var_s);
            abs_by_source(:, s, k)  = x_hat .* sqrt(var_s);
        end
        n_used = n_used + 1;
    end

    if n_used == 0
        error([newline, 'No usable pixels found during extraction (missing decomposition fields?).', newline]);
    end
    fprintf('Extracted %d usable pixels.\n', n_used);
end


function [ctr, lo, hi, idx] = quantile_bins(x, n_bins, min_count)
% Equal-count (quantile) bins along x. Returns kept-bin centers/edges and a
% per-point bin index (NaN for points not in any kept bin).
    x = x(:).';
    edges = quantile(x, linspace(0, 1, n_bins + 1));
    edges = unique(edges);                       % collapse ties
    nb = numel(edges) - 1;
    idx = nan(1, numel(x));
    ctr = nan(1, nb); lo = nan(1, nb); hi = nan(1, nb);
    keep = false(1, nb);
    for b = 1:nb
        if b < nb
            in = x >= edges(b) & x < edges(b+1);
        else
            in = x >= edges(b) & x <= edges(b+1);   % include right edge in last bin
        end
        if nnz(in) >= min_count
            idx(in) = b;
            ctr(b)  = median(x(in), 'omitnan');
            lo(b)   = edges(b);
            hi(b)   = edges(b+1);
            keep(b) = true;
        end
    end
    % compress to kept bins and remap idx to 1..n_kept
    kept_ids = find(keep);
    remap = nan(1, nb); remap(kept_ids) = 1:numel(kept_ids);
    good = ~isnan(idx);
    idx(good) = remap(idx(good));
    ctr = ctr(kept_ids); lo = lo(kept_ids); hi = hi(kept_ids);
end


function [m, plo, phi] = bin_stats(y, bin_idx, n_bins)
% Median and 16/84th percentile of y within each (kept) bin.
    m = nan(1, n_bins); plo = nan(1, n_bins); phi = nan(1, n_bins);
    for b = 1:n_bins
        vals = y(bin_idx == b);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        m(b)   = median(vals);
        pr     = prctile(vals, [16 84]);
        plo(b) = pr(1); phi(b) = pr(2);
    end
end


function shade_band(x, lo, hi, color)
% Shaded 16-84 band; skips segments where percentiles are NaN.
    good = isfinite(x) & isfinite(lo) & isfinite(hi);
    if nnz(good) < 2, return; end
    x = x(good); lo = lo(good); hi = hi(good);
    fill([x, fliplr(x)], [lo, fliplr(hi)], color, 'FaceAlpha', 0.12, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end


function maybe_save(do_save, save_directory, basename)
    if ~do_save, return; end
    fig_dir = fullfile(save_directory, 'saved_figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    stamp = char(datetime('now', 'Format', 'yyyy-MM-dd'));
    fname = fullfile(fig_dir, [basename, '_', stamp]);
    exportgraphics(gcf, [fname, '.png'], 'Resolution', 200);
    savefig(gcf, [fname, '.fig']);
    fprintf('Saved figure: %s.(png|fig)\n', fname);
end
