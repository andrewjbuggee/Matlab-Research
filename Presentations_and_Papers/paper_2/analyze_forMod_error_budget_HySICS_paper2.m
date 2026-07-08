%% Analyze and plot the forward-model-parameter error budget for the HySICS synthetic retrievals
%
% Companion to backfill_forMod_jacobian_HySICS_paper2.m, which must be run FIRST
% (on Alpine) so every retrieval carries the decomposition fields
%   forward_model_parameter_cov, forMod_sensitivity_log,
%   forMod_param_variance_log,   forMod_param_names.
% This script reads those, aggregates across the 69 synthetic VOCALS-REx pixels,
% and views the error budget five ways.
%
% CONTEXT / STRUCTURAL FACTS
% --------------------------
% The retrieval is in LOG space, so sqrt(diag(posterior_cov_log)) is the
% FRACTIONAL (relative) 1-sigma uncertainty of each retrieved variable; the
% "delta method" (x = exp(p) -> sigma_x = x_hat * sigma_p) turns it into an
% absolute uncertainty. The forward-model error was folded into the effective
% measurement covariance during the retrieval
%       S_e = S_y + K_b * S_b * K_b'
% so posterior_cov_log ALREADY contains it, and
%       S_f = G * K_b * S_b * K_b' * G'
% is a genuine sub-component of the total posterior covariance — which is what
% makes the "fraction of total" ratios below meaningful.
%
% Because S_b is DIAGONAL, S_f separates into an independent sum over the two
% physical HySICS forward-model sources:
%       r_e profile (adiabatic assumption)  and  cloud top height.
% (Unlike the EMIT retrievals, effective variance is NOT a source here.)
%
% THE FIVE VIEWS
%   FIG 1  Error-budget bar chart — median fractional 1-sigma from each source,
%          16-84th whiskers, FM-total and posterior-total overlaid.
%          -> "How big is each source, in relative terms?"
%   FIG 2  Share of VARIANCE (the statistically correct "fraction of total").
%          (a) each source as a share of the forward-model variance;
%          (b) forward-model total as a share of the total posterior variance.
%          Shares use sigma^2 because errors add in quadrature; doing it with
%          sigma directly overcounts.
%          -> "What fraction of the retrieval error is forward model, and which
%              source dominates it?"
%   FIG 3  Uncertainty vs retrieved tau_c (quantile-binned), fractional + absolute.
%          -> "Does the forward-model error grow or shrink with optical depth?"
%   FIG 4  Forward-model share of total variance vs retrieved tau_c.
%          -> "At what tau_c does the forward model start to matter?"
%   FIG 5  [SYNTHETIC-ONLY] Reliability: actual |retrieved - true| vs the
%          predicted 1-sigma, plus the standardized-error spread. Because these
%          are synthetic retrievals of known VOCALS-REx profiles, we can check
%          whether the predicted uncertainty (with the forward-model term)
%          actually covers the true error. EMIT cannot do this.
%
% tau_c is binned by QUANTILE (equal count), not fixed width, because retrieved
% tau_c is skewed; equal-count bins keep each bin's median equally determined.
%
% By Andrew John Buggee
%%

clear variables

% =============================== configuration ===============================
which_computer = whatComputer();

switch which_computer
    case 'andrewbuggee'
        retrieval_directory = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/', ...
            'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
        save_directory = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/';
    case 'anbu8374'
        retrieval_directory = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/', ...
            'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
        save_directory = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/';
    case 'curc'
        retrieval_directory = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/', ...
            'Droplet_profile_retrievals/paper2_variableSweep/', ...
            'full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
        save_directory = '/projects/anbu8374/Matlab-Research/Presentations_and_Papers/paper_2/';
    otherwise
        error([newline, 'Define retrieval_directory / save_directory for "', which_computer, '".', newline])
end

n_tau_bins    = 6;      % number of equal-count (quantile) bins along retrieved tau_c
min_bin_count = 5;      % drop tau bins with fewer than this many pixels
save_figures  = false;  % true -> write PNG/FIG into save_directory/saved_figures
save_budget   = false;  % true -> write HySICS_forMod_error_budget_<date>.mat
% =============================================================================


% ---- variable and source bookkeeping (row order of retrieval / posterior_cov_log) ----
% ln(r_top), ln(r_bot), ln(tau_c), ln(acpw)
var_fields = {'r_top', 'r_bot', 'tau_c', 'acpw'};
var_labels = {'r_{top}', 'r_{bot}', '\tau_c', 'acpw'};
var_units  = {'\mum', '\mum', '', 'cm'};
n_vars     = numel(var_fields);
i_tau      = find(strcmp(var_fields, 'tau_c'));

% two HySICS forward-model sources (must match forMod_param_names strings)
source_fields = {'re_profile', 'cloud_top_height'};
source_labels = {'r_e profile (adiabatic)', 'cloud top height'};
n_src         = numel(source_fields);

src_colors   = [0.20 0.45 0.70;    % blue   - r_e profile
                0.90 0.60 0.10];   % orange - cloud top height
col_fm_total = [0 0 0];
col_tot_post = [0.85 0.15 0.15];


%% ---- extract per-pixel arrays from the backfilled retrievals ----
[frac_total, frac_fm_total, frac_by_source, abs_total, abs_fm_total, abs_by_source, ...
 xhat_all, xtrue_all, n_used, n_missing] = ...
    extract_hysics_arrays(retrieval_directory, source_fields, var_fields, n_vars, n_src);

if n_missing > 0
    fprintf(['\n%d file(s) lacked the decomposition fields. Run ', ...
        'backfill_forMod_jacobian_HySICS_paper2.m first (on Alpine).\n'], n_missing);
end

% keep only fully finite pixels for the budget
valid = all(isfinite(xhat_all), 1) & all(isfinite(frac_total), 1) & all(isfinite(frac_fm_total), 1);
frac_total = frac_total(:, valid);        frac_fm_total = frac_fm_total(:, valid);
frac_by_source = frac_by_source(:, :, valid);
abs_total = abs_total(:, valid);          abs_fm_total = abs_fm_total(:, valid);
abs_by_source = abs_by_source(:, :, valid);
xhat_all = xhat_all(:, valid);            xtrue_all = xtrue_all(:, valid);
N = size(xhat_all, 2);

if N == 0
    error([newline, 'No usable pixels with decomposition fields found. Backfill first.', newline])
end
fprintf('\nAnalyzing %d valid synthetic pixels.\n', N);
if N < min_bin_count
    warning('Only %d valid pixels; tau-binned figures (3,4) may be unreliable.', N);
end


%% ---- derived quantities, all in VARIANCE (sigma^2) space ----
var_total     = frac_total.^2;                                  % n_vars x N
var_fm_total  = frac_fm_total.^2;                               % n_vars x N
var_by_source = frac_by_source.^2;                              % n_vars x n_src x N

share_src_of_fm   = var_by_source ./ max(sum(var_by_source, 2), eps);   % source share of FM budget
share_fm_of_total = var_fm_total ./ max(var_total, eps);                % FM share of posterior

tau_c_retrieved = xhat_all(i_tau, :);                           % 1 x N


%% =====================================================================
%  FIGURE 1 — error-budget bar chart (fractional 1-sigma, %)
%  =====================================================================
med_src_frac = 100 * median(frac_by_source, 3, 'omitnan');      % n_vars x n_src
p_src_frac   = 100 * prctile(frac_by_source, [16 84], 3);       % n_vars x n_src x 2
med_fm_frac  = 100 * median(frac_fm_total, 2, 'omitnan');       % n_vars x 1
med_tot_frac = 100 * median(frac_total,    2, 'omitnan');       % n_vars x 1

figure('Position', [80 80 900 520]);
b = bar(med_src_frac, 'grouped');
for s = 1:n_src
    b(s).FaceColor = src_colors(s, :);
    b(s).DisplayName = source_labels{s};
end
hold on
for s = 1:n_src
    xc = b(s).XEndPoints;
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
title(sprintf('HySICS synthetic forward-model error budget (median of %d pixels; whiskers 16-84th pct)', N));
legend([b, h_fm, h_tot], 'Location', 'best', 'Interpreter', 'tex');
grid on; grid minor
maybe_save(save_figures, save_directory, 'HySICS_forMod_errorBudget_fig1_bar');


%% =====================================================================
%  FIGURE 2 — share of variance (the correct "fraction of total")
%  =====================================================================
figure('Position', [110 110 1000 460]);

subplot(1, 2, 1);
med_share_src_of_fm = 100 * median(share_src_of_fm, 3, 'omitnan');
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

subplot(1, 2, 2);
med_share_fm_of_total = 100 * median(share_fm_of_total, 2, 'omitnan');
p_share_fm_of_total   = 100 * prctile(share_fm_of_total, [16 84], 2);
bar(1:n_vars, med_share_fm_of_total, 0.6, 'FaceColor', [0.5 0.5 0.5]);
hold on
errorbar(1:n_vars, med_share_fm_of_total, ...
    med_share_fm_of_total - p_share_fm_of_total(:, 1), p_share_fm_of_total(:, 2) - med_share_fm_of_total, ...
    'LineStyle', 'none', 'Color', 'k', 'CapSize', 6);
set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
ylabel('Forward-model share of total variance (%)', 'Interpreter', 'tex');
title('(b) How much of the total retrieval error is forward model');
grid on
maybe_save(save_figures, save_directory, 'HySICS_forMod_errorBudget_fig2_shares');


%% =====================================================================
%  FIGURE 3 — uncertainty vs retrieved tau_c (quantile-binned)
%  top row = fractional (%), bottom row = absolute (retrieved units)
%  =====================================================================
[bin_ctr, ~, ~, bin_idx] = quantile_bins(tau_c_retrieved, n_tau_bins, min_bin_count);
n_bins_kept = numel(bin_ctr);

figure('Position', [140 60 1200 720]);
for i = 1:n_vars
    subplot(2, n_vars, i); hold on
    lg = gobjects(0);
    for s = 1:n_src
        y = squeeze(frac_by_source(i, s, :)).' * 100;
        [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
        shade_band(bin_ctr, plo, phi, src_colors(s, :));
        lg(end+1) = plot(bin_ctr, m, '-o', 'Color', src_colors(s, :), 'MarkerFaceColor', src_colors(s, :), ...
            'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', source_labels{s}); %#ok<SAGROW>
    end
    m_fm  = bin_stats(frac_fm_total(i, :) * 100, bin_idx, n_bins_kept);
    m_tot = bin_stats(frac_total(i, :)    * 100, bin_idx, n_bins_kept);
    lg(end+1) = plot(bin_ctr, m_fm,  's--', 'Color', col_fm_total, 'MarkerFaceColor', col_fm_total, ...
        'MarkerSize', 4, 'DisplayName', 'forward-model (all)'); %#ok<SAGROW>
    lg(end+1) = plot(bin_ctr, m_tot, 'd--', 'Color', col_tot_post, 'MarkerFaceColor', col_tot_post, ...
        'MarkerSize', 4, 'DisplayName', 'total posterior'); %#ok<SAGROW>
    xlabel('retrieved \tau_c', 'Interpreter', 'tex');
    ylabel('fractional 1\sigma (%)', 'Interpreter', 'tex');
    title(var_labels{i}, 'Interpreter', 'tex'); grid on
    if i == 1, legend(lg, 'Location', 'best', 'Interpreter', 'tex', 'FontSize', 7); end

    subplot(2, n_vars, n_vars + i); hold on
    for s = 1:n_src
        y = squeeze(abs_by_source(i, s, :)).';
        [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
        shade_band(bin_ctr, plo, phi, src_colors(s, :));
        plot(bin_ctr, m, '-o', 'Color', src_colors(s, :), 'MarkerFaceColor', src_colors(s, :), ...
            'MarkerSize', 4, 'LineWidth', 1.5);
    end
    m_fm_a  = bin_stats(abs_fm_total(i, :), bin_idx, n_bins_kept);
    m_tot_a = bin_stats(abs_total(i, :),    bin_idx, n_bins_kept);
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
maybe_save(save_figures, save_directory, 'HySICS_forMod_errorBudget_fig3_vsTau');


%% =====================================================================
%  FIGURE 4 — forward-model share of total variance vs retrieved tau_c
%  =====================================================================
figure('Position', [170 90 780 520]); hold on
line_cols = lines(n_vars);
for i = 1:n_vars
    y = share_fm_of_total(i, :) * 100;
    [m, plo, phi] = bin_stats(y, bin_idx, n_bins_kept);
    shade_band(bin_ctr, plo, phi, line_cols(i, :));
    plot(bin_ctr, m, '-o', 'Color', line_cols(i, :), 'MarkerFaceColor', line_cols(i, :), ...
        'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName', var_labels{i});
end
xlabel('retrieved \tau_c', 'Interpreter', 'tex');
ylabel('forward-model share of total variance (%)', 'Interpreter', 'tex');
title(sprintf('When does the forward model matter? (%d quantile bins; shaded 16-84th pct)', n_bins_kept));
legend('Location', 'best', 'Interpreter', 'tex'); grid on; grid minor
maybe_save(save_figures, save_directory, 'HySICS_forMod_errorBudget_fig4_shareVsTau');


%% =====================================================================
%  FIGURE 5 — [SYNTHETIC-ONLY] reliability: actual error vs predicted 1-sigma
%  =====================================================================
have_truth = any(isfinite(xtrue_all(:)));
if have_truth
    actual_err = abs(xhat_all - xtrue_all);      % n_vars x N, absolute
    stderr     = (xhat_all - xtrue_all) ./ max(abs_total, eps);   % standardized (error / predicted sigma)

    figure('Position', [200 60 1200 520]);
    for i = 1:n_vars
        subplot(2, n_vars, i); hold on
        ok = isfinite(actual_err(i, :)) & isfinite(abs_total(i, :));
        if nnz(ok) >= 1
            % predicted 1-sigma WITH the forward-model term (abs_total) and the
            % measurement-only part (abs_total minus FM in quadrature) for contrast
            abs_meas_only = sqrt(max(abs_total(i, :).^2 - abs_fm_total(i, :).^2, 0));
            scatter(abs_total(i, ok), actual_err(i, ok), 22, [0.2 0.4 0.7], 'filled', ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', 'predicted (with FM)');
            scatter(abs_meas_only(ok), actual_err(i, ok), 16, [0.7 0.7 0.7], ...
                'MarkerEdgeColor', [0.5 0.5 0.5], 'DisplayName', 'predicted (meas only)');
            mx = max([abs_total(i, ok), actual_err(i, ok)], [], 'omitnan');
            plot([0 mx], [0 mx], 'k--', 'DisplayName', '1:1');
        end
        xlabel(sprintf('predicted 1\\sigma%s', ternary(isempty(var_units{i}), '', [' (' var_units{i} ')'])), ...
            'Interpreter', 'tex');
        ylabel('actual |retrieved-true|', 'Interpreter', 'tex');
        title(var_labels{i}, 'Interpreter', 'tex'); grid on; axis equal
        if i == 1, legend('Location', 'best', 'FontSize', 7); end

        % standardized-error spread: well-calibrated -> ~ N(0,1), |z|<=1 ~68%
        subplot(2, n_vars, n_vars + i); hold on
        z = stderr(i, isfinite(stderr(i, :)));
        if ~isempty(z)
            histogram(z, 'Normalization', 'pdf', 'FaceColor', [0.2 0.4 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            xx = linspace(-4, 4, 200);
            plot(xx, exp(-xx.^2/2)/sqrt(2*pi), 'k-', 'LineWidth', 1.2);
            frac_in = mean(abs(z) <= 1);
            title(sprintf('%s: %.0f%% within \\pm1\\sigma', var_labels{i}, 100*frac_in), 'Interpreter', 'tex');
        end
        xlabel('(retrieved-true)/predicted 1\sigma', 'Interpreter', 'tex');
        ylabel('pdf'); xlim([-4 4]); grid on
    end
    sgtitle(sprintf('Reliability of the predicted uncertainty (%d synthetic pixels; top: does it cover the true error?)', N));
    maybe_save(save_figures, save_directory, 'HySICS_forMod_errorBudget_fig5_reliability');
else
    fprintf('\n[FIG 5 skipped] No truth values recovered from GN_inputs.measurement.\n');
end


%% ---- concise text summary ----
fprintf('\n==================== HySICS synthetic forward-model error budget (median of %d pixels) ====================\n', N);
for i = 1:n_vars
    fprintf('\n%-6s (retrieved %s median = %.3g)\n', var_fields{i}, var_labels{i}, median(xhat_all(i, :), 'omitnan'));
    fprintf('   forward-model 1sigma : %6.2f%% of value   |   %6.1f%% of total variance\n', ...
        med_fm_frac(i), 100*median(share_fm_of_total(i, :), 'omitnan'));
    for s = 1:n_src
        fprintf('     - %-22s : %6.2f%% of value   |   %6.1f%% of FM variance\n', ...
            source_labels{s}, med_src_frac(i, s), 100*median(share_src_of_fm(i, s, :), 'omitnan'));
    end
end
fprintf('\n(Shares are variance-based; the two sources sum to ~100%% of the forward-model variance.)\n\n');


%% ---- optional: save the aggregated budget ----
if save_budget
    budget.n_pixels      = N;
    budget.var_fields    = var_fields;
    budget.source_fields = source_fields;
    budget.frac_total    = frac_total;
    budget.frac_fm_total = frac_fm_total;
    budget.frac_by_source = frac_by_source;
    budget.abs_total     = abs_total;
    budget.abs_fm_total  = abs_fm_total;
    budget.abs_by_source = abs_by_source;
    budget.xhat_all      = xhat_all;
    budget.xtrue_all     = xtrue_all;
    out_name = fullfile(save_directory, ['HySICS_forMod_error_budget_', char(datetime('today')), '.mat']);
    save(out_name, 'budget', 'retrieval_directory');
    fprintf('Saved budget to:\n  %s\n', out_name);
end


%% ============================ local functions ============================

function [frac_total, frac_fm_total, frac_by_source, abs_total, abs_fm_total, abs_by_source, ...
          xhat_all, xtrue_all, n_used, n_missing] = ...
          extract_hysics_arrays(retrieval_directory, source_fields, var_fields, n_vars, n_src)
% Read every backfilled HySICS retrieval and build the per-pixel arrays.

    file_list = dir(fullfile(retrieval_directory, 'dropletRetrieval_HySICS_*.mat'));
    if isempty(file_list)
        error([newline, 'No dropletRetrieval_HySICS_*.mat files in:', newline, retrieval_directory, newline]);
    end
    n_files = numel(file_list);

    frac_total     = nan(n_vars, n_files);
    frac_fm_total  = nan(n_vars, n_files);
    frac_by_source = nan(n_vars, n_src, n_files);
    abs_total      = nan(n_vars, n_files);
    abs_fm_total   = nan(n_vars, n_files);
    abs_by_source  = nan(n_vars, n_src, n_files);
    xhat_all       = nan(n_vars, n_files);
    xtrue_all      = nan(n_vars, n_files);
    n_used = 0; n_missing = 0;

    needed = {'retrieval', 'posterior_cov_log', 'forward_model_parameter_cov', ...
        'forMod_sensitivity_log', 'forMod_param_variance_log', 'forMod_param_names'};

    for k = 1:n_files
        S = load(fullfile(file_list(k).folder, file_list(k).name), 'GN_outputs', 'GN_inputs');
        if ~isfield(S, 'GN_outputs') || isempty(S.GN_outputs), n_missing = n_missing + 1; continue; end
        GN = S.GN_outputs;
        if ~all(isfield(GN, needed)), n_missing = n_missing + 1; continue; end

        x_hat = GN.retrieval(:, end);
        if numel(x_hat) ~= n_vars || any(~isfinite(x_hat)), continue; end

        % total posterior: fractional (log-space std) and absolute (delta method)
        sig_frac_total = sqrt(diag(GN.posterior_cov_log));
        frac_total(:, k) = sig_frac_total;
        abs_total(:, k)  = x_hat .* sig_frac_total;
        xhat_all(:, k)   = x_hat;

        % forward-model total
        sig_frac_fm = sqrt(diag(GN.forward_model_parameter_cov));
        frac_fm_total(:, k) = sig_frac_fm;
        abs_fm_total(:, k)  = x_hat .* sig_frac_fm;

        % per-source decomposition: Var(ln x_i) from source s = sum_j Var(ln b_j)*sens_ij^2
        sens  = GN.forMod_sensitivity_log;              % n_vars x p
        var_b = GN.forMod_param_variance_log(:).';      % 1 x p
        names = GN.forMod_param_names;                  % 1 x p
        contrib_var = (sens.^2) .* var_b;               % n_vars x p
        for s = 1:n_src
            cols  = strcmp(names, source_fields{s});
            var_s = sum(contrib_var(:, cols), 2);
            frac_by_source(:, s, k) = sqrt(var_s);
            abs_by_source(:, s, k)  = x_hat .* sqrt(var_s);
        end

        % truth (synthetic advantage) from GN_inputs.measurement, best-effort
        if isfield(S, 'GN_inputs') && isfield(S.GN_inputs, 'measurement')
            xtrue_all(:, k) = extract_truth(S.GN_inputs.measurement, var_fields);
        end

        n_used = n_used + 1;
    end

    fprintf('Extracted %d usable pixels (%d missing decomposition/empty).\n', n_used, n_missing);
end


function xt = extract_truth(meas, var_fields)
% Recover the true [r_top; r_bot; tau_c; acpw] used to generate the synthetic
% measurement. tau_c and acpw are stored directly; r_top/r_bot are the cloud-top
% and cloud-base values of the in-situ effective-radius profile. Missing -> NaN.
    xt = nan(numel(var_fields), 1);
    if isfield(meas, 're_prof') && isfield(meas, 'z') && ~isempty(meas.re_prof)
        re = meas.re_prof(:); z = meas.z(:);
        if numel(re) == numel(z) && numel(z) >= 2
            [~, itop] = max(z); [~, ibot] = min(z);
            xt(strcmp(var_fields, 'r_top')) = re(itop);
            xt(strcmp(var_fields, 'r_bot')) = re(ibot);
        end
    end
    if isfield(meas, 'tau_c'), xt(strcmp(var_fields, 'tau_c')) = meas.tau_c; end
    if isfield(meas, 'actpw'), xt(strcmp(var_fields, 'acpw'))  = meas.actpw; end
end


function [ctr, lo, hi, idx] = quantile_bins(x, n_bins, min_count)
% Equal-count (quantile) bins along x. Returns kept-bin centers/edges and a
% per-point bin index (NaN for points not in any kept bin).
    x = x(:).';
    edges = quantile(x, linspace(0, 1, n_bins + 1));
    edges = unique(edges);
    nb = numel(edges) - 1;
    idx = nan(1, numel(x));
    ctr = nan(1, nb); lo = nan(1, nb); hi = nan(1, nb);
    keep = false(1, nb);
    for b = 1:nb
        if b < nb
            in = x >= edges(b) & x < edges(b+1);
        else
            in = x >= edges(b) & x <= edges(b+1);
        end
        if nnz(in) >= min_count
            idx(in) = b; ctr(b) = median(x(in), 'omitnan');
            lo(b) = edges(b); hi(b) = edges(b+1); keep(b) = true;
        end
    end
    kept_ids = find(keep);
    remap = nan(1, nb); remap(kept_ids) = 1:numel(kept_ids);
    good = ~isnan(idx); idx(good) = remap(idx(good));
    ctr = ctr(kept_ids); lo = lo(kept_ids); hi = hi(kept_ids);
end


function [m, plo, phi] = bin_stats(y, bin_idx, n_bins)
% Median and 16/84th percentile of y within each kept bin.
    m = nan(1, n_bins); plo = nan(1, n_bins); phi = nan(1, n_bins);
    for b = 1:n_bins
        vals = y(bin_idx == b); vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        m(b) = median(vals);
        pr = prctile(vals, [16 84]); plo(b) = pr(1); phi(b) = pr(2);
    end
end


function shade_band(x, lo, hi, color)
% Shaded 16-84 band; skips NaN segments.
    good = isfinite(x) & isfinite(lo) & isfinite(hi);
    if nnz(good) < 2, return; end
    x = x(good); lo = lo(good); hi = hi(good);
    fill([x, fliplr(x)], [lo, fliplr(hi)], color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
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
