%% Analyze and plot the forward-model-parameter error budget for the EMIT-Aqua retrievals
%
% Visualization / analysis companion to compute_forMod_error_budget_paper2.m.
% For every retrieved pixel the forward-model (FM) parameter uncertainty has
% been decomposed, in log space, into three physical sources:
%       r_e profile (adiabatic assumption), cloud top height, effective variance.
% This script views that decomposition two ways plus a relative-magnitude view:
%
%   FIG 1  Fraction of the TOTAL retrieval uncertainty from each FM source
%          Per retrieved variable, a 100%-stacked bar splitting the total
%          posterior VARIANCE into the three FM sources plus the remainder
%          (measurement noise + a priori). This is the direct answer to
%          "what fraction of the total uncertainty is due to each FM assumption?"
%
%   FIG 2  Absolute 1-sigma contribution of each FM source (physical units)
%          Per retrieved variable, grouped bars of each source's 1-sigma in the
%          variable's own units (r_top, r_bot in microns; tau_c dimensionless;
%          acpw in cm), with the FM-total and posterior-total overlaid.
%
%   FIG 3  Relative (fractional) 1-sigma error budget
%          Same decomposition expressed as a percent of the retrieved value.
%
% WHY VARIANCE, NOT SIGMA, FOR "FRACTION OF TOTAL":
%   The retrieval folds the FM error into the effective measurement covariance
%   (calc_retrieval_..._perPixel.m ~line 1101):
%       total_meas_fm_cov = measurement_cov_log + K_b * S_b * K_b'
%   so the total posterior variance is an ADDITIVE budget
%       sigma_total^2 = sigma_meas+prior^2 + sigma_fm^2 ,   sigma_fm^2 = sum_s sigma_s^2
%   (S_b is diagonal, so the FM sources add in variance). Fractions of the total
%   are therefore computed as sigma_s^2 / sigma_total^2. Doing it with sigma
%   directly would overcount because the pieces do not add linearly.
%   Absolute 1-sigma values in FIG 2 also combine in quadrature, which is why
%   they are shown as GROUPED (not stacked) bars.
%
% Data source: prefers the newest EMIT_Aqua_forMod_error_budget_*.mat saved by
% compute_forMod_error_budget_paper2.m; if none is found it runs the identical
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

save_figures = false;    % set true to write PNG/FIG into save_directory/saved_figures
% ----------------------------------------------------------------------


% ---- variable and source bookkeeping (must match the compute script) ----
% Row order of the retrieval / posterior_cov_log: ln(r_top), ln(r_bot), ln(tau_c), ln(acpw)
var_fields = {'r_top', 'r_bot', 'tau_c', 'acpw'};
var_labels = {'r_{top}', 'r_{bot}', '\tau_c', 'acpw'};
var_units  = {'\mum', '\mum', '', 'cm'};      % absolute-unit labels for retrieved values
n_vars     = numel(var_fields);

source_fields = {'re_profile', 'cloud_top_height', 'effective_variance'};
source_labels = {'r_e profile (adiabatic)', 'cloud top height', 'effective variance'};
n_src         = numel(source_fields);

% consistent colors: 3 FM sources, a gray for the measurement+prior remainder,
% then FM-total (black) and posterior-total (red) as overlaid markers
src_colors = [0.20 0.45 0.70;    % blue   - r_e profile
              0.90 0.60 0.10;    % orange - cloud top height
              0.35 0.70 0.35];   % green  - effective variance
col_remain   = [0.75 0.75 0.75]; % gray   - measurement noise + a priori
col_fm_total = [0 0 0];
col_tot_post = [0.85 0.15 0.15];


%% ---- get the per-pixel arrays (load saved, else extract) ----
[frac_total, frac_fm_total, frac_by_source, ...
 abs_total,  abs_fm_total,  abs_by_source, xhat_all, n_used, src_from_data] = ...
    get_perpixel_arrays(save_directory, retrieval_directory, source_fields, n_vars, n_src);

if ~isempty(src_from_data) && ~isequal(src_from_data(:)', source_fields)
    warning(['Source-field order in the data (%s) differs from this script''s ', ...
        'assumption (%s). Check source_fields.'], strjoin(src_from_data, ', '), strjoin(source_fields, ', '));
end

% ---- keep only fully finite pixels ----
src_finite = reshape(all(all(isfinite(frac_by_source), 1), 2), 1, []);   % 1 x N
valid = all(isfinite(xhat_all), 1) & all(isfinite(frac_total), 1) & all(isfinite(frac_fm_total), 1) ...
      & src_finite;
frac_total     = frac_total(:, valid);
frac_fm_total  = frac_fm_total(:, valid);
frac_by_source = frac_by_source(:, :, valid);
abs_total      = abs_total(:, valid);
abs_fm_total   = abs_fm_total(:, valid);
abs_by_source  = abs_by_source(:, :, valid);
xhat_all       = xhat_all(:, valid);
N              = size(xhat_all, 2);

fprintf('\nAnalyzing %d valid pixels (of %d used in extraction).\n', N, n_used);


%% ---- derived per-pixel quantities (VARIANCE space is the additive one) ----
var_total     = frac_total.^2;                 % n_vars x N        total posterior
var_fm_total  = frac_fm_total.^2;              % n_vars x N        all FM sources
var_by_source = frac_by_source.^2;             % n_vars x n_src x N per source

% remainder of the total NOT explained by the forward model (measurement + prior).
% max(.,0) guards against tiny negative round-off; verified >=0 in the data.
var_remaining = max(var_total - var_fm_total, 0);              % n_vars x N

% ---- fraction of the TOTAL variance carried by each component (per pixel) ----
share_src_of_total = var_by_source ./ reshape(max(var_total, eps), n_vars, 1, N);  % n_vars x n_src x N
share_remaining    = var_remaining ./ max(var_total, eps);                          % n_vars x N
share_fm_of_total  = var_fm_total  ./ max(var_total, eps);                          % n_vars x N

% flag any pixel where the FM variance exceeds the total (should be none)
n_over = nnz(share_fm_of_total > 1 + 1e-6);
if n_over > 0
    warning('%d pixel-variables have FM variance > total variance (check the covariance math).', n_over);
end


%% =====================================================================
%  FIGURE 1 — fraction of the TOTAL retrieval variance from each source
%  (100%-stacked; the direct answer to "what fraction is each FM assumption")
%  =====================================================================
% median share of each component across pixels, then renormalize the four
% medians to sum to exactly 100% for a clean stacked bar (medians are not
% additive; the renormalization is cosmetic and <1% in practice).
med_share_src    = median(share_src_of_total, 3, 'omitnan');   % n_vars x n_src
med_share_remain = median(share_remaining, 2, 'omitnan');      % n_vars x 1
stack_raw  = [med_share_src, med_share_remain];                % n_vars x (n_src+1)
stack_pct  = 100 * stack_raw ./ sum(stack_raw, 2);             % renormalized to 100%

figure('Position', [80 80 900 540]);
bS = bar(stack_pct, 'stacked');
stack_colors = [src_colors; col_remain];
stack_names  = [source_labels, {'measurement + a priori'}];
for c = 1:numel(bS)
    bS(c).FaceColor = stack_colors(c, :);
    bS(c).DisplayName = stack_names{c};
end
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Share of total retrieval variance (%)', 'Interpreter', 'tex');
ylim([0 100]);
title(sprintf('Fraction of total retrieval uncertainty by source (median of %d pixels)', N));
legend(bS, 'Location', 'eastoutside', 'Interpreter', 'tex');
grid on
% annotate the forward-model total share (sum of the three source segments)
for i = 1:n_vars
    fm_pct = sum(stack_pct(i, 1:n_src));
    text(i, min(fm_pct + 3, 98), sprintf('%.0f%% FM', fm_pct), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig1_fractionOfTotal');


%% =====================================================================
%  FIGURE 2 — absolute 1-sigma contribution in physical units (grouped)
%  one panel per retrieved variable, because the units differ
%  =====================================================================
med_abs_src  = median(abs_by_source, 3, 'omitnan');           % n_vars x n_src
p_abs_src    = prctile(abs_by_source, [16 84], 3);            % n_vars x n_src x 2
med_abs_fm   = median(abs_fm_total, 2, 'omitnan');            % n_vars x 1
med_abs_tot  = median(abs_total,    2, 'omitnan');            % n_vars x 1

figure('Position', [110 90 1050 640]);
for i = 1:n_vars
    subplot(2, 2, i); hold on
    hb = bar(1:n_src, med_abs_src(i, :), 0.65, 'FaceColor', 'flat');
    hb.CData = src_colors;
    % 16-84 whiskers
    lo = squeeze(p_abs_src(i, :, 1)).';   hi = squeeze(p_abs_src(i, :, 2)).';
    errorbar(1:n_src, med_abs_src(i, :), med_abs_src(i, :) - lo.', hi.' - med_abs_src(i, :), ...
        'LineStyle', 'none', 'Color', 0.25*[1 1 1], 'CapSize', 6);
    % FM-total and posterior-total as horizontal reference lines (labels on the
    % right so they clear the top-left panel label)
    yline(med_abs_fm(i),  '--', 'FM total', 'Color', col_fm_total, ...
        'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'Interpreter', 'tex');
    yline(med_abs_tot(i), '-',  'posterior total', 'Color', col_tot_post, 'LineWidth', 1.2, ...
        'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'Interpreter', 'tex');
    ylim([0, 1.22 * med_abs_tot(i)]);   % headroom for the in-panel label
    set(gca, 'XTick', 1:n_src, 'XTickLabel', {'r_e prof', 'CTH', 'eff var'});
    if isempty(var_units{i})
        ylabel('absolute 1\sigma', 'Interpreter', 'tex');
        ttl = sprintf('%s  (median = %.3g)', var_labels{i}, median(xhat_all(i, :), 'omitnan'));
    else
        ylabel(sprintf('absolute 1\\sigma (%s)', var_units{i}), 'Interpreter', 'tex');
        ttl = sprintf('%s  (median = %.3g %s)', var_labels{i}, median(xhat_all(i, :), 'omitnan'), var_units{i});
    end
    title(ttl, 'Interpreter', 'tex');
    grid on
end
sgtitle(sprintf(['Absolute 1\\sigma forward-model contribution (median of %d pixels; whiskers 16-84th pct)', ...
    '   —   sources add in quadrature: \\sigma_{FM}^2 = \\Sigma_s \\sigma_s^2'], N), 'Interpreter', 'tex', 'FontSize', 11);
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig2_absoluteUnits');


%% =====================================================================
%  FIGURE 3 — relative (fractional) 1-sigma error budget (grouped, %)
%  =====================================================================
med_src_frac = 100 * median(frac_by_source, 3, 'omitnan');    % n_vars x n_src (%)
p_src_frac   = 100 * prctile(frac_by_source, [16 84], 3);     % n_vars x n_src x 2 (%)
med_fm_frac  = 100 * median(frac_fm_total, 2, 'omitnan');     % n_vars x 1
med_tot_frac = 100 * median(frac_total,    2, 'omitnan');     % n_vars x 1

figure('Position', [140 100 900 540]);
b = bar(med_src_frac, 'grouped');
for s = 1:n_src
    b(s).FaceColor = src_colors(s, :);
    b(s).DisplayName = source_labels{s};
end
hold on
for s = 1:n_src
    xc = b(s).XEndPoints;                              % bar x-centers for this source
    lo = p_src_frac(:, s, 1).';                        % already in % (single scaling)
    hi = p_src_frac(:, s, 2).';
    errorbar(xc, med_src_frac(:, s).', med_src_frac(:, s).' - lo, hi - med_src_frac(:, s).', ...
        'LineStyle', 'none', 'Color', 0.3*[1 1 1], 'CapSize', 4, 'HandleVisibility', 'off');
end
h_fm  = plot(1:n_vars, med_fm_frac,  's', 'MarkerFaceColor', col_fm_total, 'MarkerEdgeColor', col_fm_total, ...
    'MarkerSize', 9, 'LineStyle', 'none', 'DisplayName', 'forward-model (all sources)');
h_tot = plot(1:n_vars, med_tot_frac, 'd', 'MarkerFaceColor', col_tot_post, 'MarkerEdgeColor', col_tot_post, ...
    'MarkerSize', 9, 'LineStyle', 'none', 'DisplayName', 'total posterior');
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Fractional 1\sigma uncertainty (% of retrieved value)', 'Interpreter', 'tex');
title(sprintf('Relative forward-model error budget (median of %d pixels; whiskers 16-84th pct)', N));
legend([b, h_fm, h_tot], 'Location', 'best', 'Interpreter', 'tex');
grid on
maybe_save(save_figures, save_directory, 'forMod_errorBudget_fig3_fractional');


%% ---- text summary answering the two headline questions ----
fprintf('\n================= forward-model error budget (median of %d pixels) =================\n', N);
fprintf('%-6s | %-14s | %-32s | %-30s\n', 'var', 'retrieved', 'FRACTION of total variance (%)', 'ABSOLUTE 1-sigma');
fprintf('%s\n', repmat('-', 1, 92));
for i = 1:n_vars
    u = var_units{i}; if isempty(u), u = '(-)'; end
    fprintf('%-6s | %8.3g %-5s | FM total = %5.1f%%\n', var_fields{i}, median(xhat_all(i, :), 'omitnan'), u, ...
        100*median(share_fm_of_total(i, :), 'omitnan'));
    for s = 1:n_src
        fprintf('       |                | %-20s %5.1f%% | %8.3g %s\n', source_labels{s}, ...
            100*median(share_src_of_total(i, s, :), 'omitnan'), med_abs_src(i, s), u);
    end
    fprintf('       |                | %-20s %5.1f%% | %8.3g %s\n', 'measurement + prior', ...
        100*median(share_remaining(i, :), 'omitnan'), med_abs_tot(i), u);
end
fprintf('\n(Shares are variance-based and sum to 100%% per variable; absolute values combine in quadrature.)\n\n');


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
