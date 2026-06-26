%% Forward-model-parameter error budget for the EMIT-Aqua droplet retrievals
%
% Run this AFTER all EMIT-Aqua single-pixel retrievals are saved. For every
% retrieved pixel it:
%
%   (1) transforms the log-space retrieval covariances to LINEAR space using
%       first-order error propagation (the "delta method"):
%           x = exp(p)  ->  sigma_x = x_hat * sigma_p
%       so sqrt(diag(posterior_cov_log)) is the FRACTIONAL uncertainty, and
%       multiplying by the retrieved value x_hat gives the absolute uncertainty;
%
%   (2) decomposes the forward-model-parameter contribution to the retrieval
%       uncertainty into the individual physical sources accounted for in the
%       forward model covariance matrix:
%           - r_e profile (adiabatic droplet-profile assumption)
%           - cloud top height
%           - droplet distribution effective variance
%       Because the forward model covariance S_b is DIAGONAL, the total
%       forward-model error covariance S_f = G*K_b*S_b*K_b'*G' separates into
%       an independent sum over parameters j:
%           Var(ln x_i) from param j  =  Var(ln b_j) * [ (G*K_b)_{ij} ]^2
%       The r_e-profile contribution is the sum of the per-layer variances
%       (the layers are modelled as independent in S_b).
%
%   (3) aggregates each quantity across all pixels (median + 16th/84th
%       percentiles, with the mean reported for reference).
%
% Requires the retrieval .mat files to have been produced by
% calc_retrieval_gauss_newton_EMIT_ver4_log_forMo_uncert_perPixel.m with the
% fields: posterior_cov_log, retrieval, forward_model_parameter_cov,
% forMod_sensitivity_log, forMod_param_variance_log, forMod_param_names.
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
    error([newline, 'Define retrieval_directory for this computer.', newline])
end

% central tendency used for the headline numbers: 'median' (robust, recommended)
% or 'mean'. The script computes both regardless; this only sets the figure/print.
central_stat   = 'median';

plot_figures   = true;
save_results   = true;
% ----------------------------------------------------------------------


% ---- retrieved-variable and error-source bookkeeping ----
% Row order of the retrieval / posterior_cov_log: ln(r_top), ln(r_bot), ln(tau_c), ln(acpw)
var_fields = {'r_top', 'r_bot', 'tau_c', 'acpw'};
var_labels = {'r_{top}', 'r_{bot}', '\tau_c', 'acpw'};
n_vars     = numel(var_fields);

% Physical forward-model error sources (must match the strings written by
% calc_retrieval_..._perPixel.m into GN_outputs.forMod_param_names)
source_fields = {'re_profile', 'cloud_top_height', 'effective_variance'};
source_labels = {'r_e profile (adiabatic)', 'cloud top height', 'effective variance'};
n_src         = numel(source_fields);


%% ---- find the retrieval files ----
file_list = dir(fullfile(retrieval_directory, 'EMIT_dropRetrieval_*.mat'));

if isempty(file_list)
    error([newline, 'No EMIT_dropRetrieval_*.mat files found in:', newline, ...
        retrieval_directory, newline])
end

n_files = numel(file_list);
fprintf('\nFound %d retrieval files in %s\n', n_files, retrieval_directory);


%% ---- preallocate per-pixel accumulators (NaN = "this pixel skipped") ----
% fractional (relative, dimensionless) 1-sigma uncertainties
frac_total      = nan(n_vars, n_files);            % total posterior
frac_fm_total   = nan(n_vars, n_files);            % all forward-model params combined
frac_by_source  = nan(n_vars, n_src, n_files);     % per physical source

% absolute (linear retrieval units) 1-sigma uncertainties
abs_total       = nan(n_vars, n_files);
abs_fm_total    = nan(n_vars, n_files);
abs_by_source   = nan(n_vars, n_src, n_files);

xhat_all        = nan(n_vars, n_files);            % retrieved (linear) values
decomp_residual = nan(n_vars, n_files);            % sanity-check: source sum vs S_f diag

n_skipped_empty  = 0;   % masked / empty GN_outputs
n_skipped_fields = 0;   % missing the new decomposition fields (old retrievals)
n_used           = 0;


%% ---- loop over pixels ----
for k = 1:n_files

    S = load(fullfile(file_list(k).folder, file_list(k).name), 'GN_outputs');

    if ~isfield(S, 'GN_outputs') || isempty(S.GN_outputs)
        n_skipped_empty = n_skipped_empty + 1;
        continue
    end
    GN = S.GN_outputs;

    % require all the pieces we need
    needed = {'retrieval', 'posterior_cov_log', 'forward_model_parameter_cov', ...
        'forMod_sensitivity_log', 'forMod_param_variance_log', 'forMod_param_names'};
    if ~all(isfield(GN, needed))
        n_skipped_fields = n_skipped_fields + 1;
        continue
    end

    % retrieved (linear) state vector at convergence: [r_top; r_bot; tau_c; acpw]
    x_hat = GN.retrieval(:, end);
    if numel(x_hat) ~= n_vars || any(~isfinite(x_hat))
        n_skipped_empty = n_skipped_empty + 1;
        continue
    end

    % ---- (1) total retrieval uncertainty ----
    % log-space std = FRACTIONAL uncertainty; * x_hat = absolute (delta method)
    sigma_frac_total = sqrt(diag(GN.posterior_cov_log));        % n_vars x 1
    frac_total(:, k) = sigma_frac_total;
    abs_total(:, k)  = x_hat .* sigma_frac_total;
    xhat_all(:, k)   = x_hat;

    % ---- forward-model total (all parameters combined) ----
    sigma_frac_fm = sqrt(diag(GN.forward_model_parameter_cov));  % n_vars x 1
    frac_fm_total(:, k) = sigma_frac_fm;
    abs_fm_total(:, k)  = x_hat .* sigma_frac_fm;

    % ---- (2) per-parameter -> per-source decomposition ----
    % variance contribution of parameter j to ln(x_i):  Var(ln b_j) * sens_ij^2
    sens   = GN.forMod_sensitivity_log;                  % n_vars x p
    var_b  = GN.forMod_param_variance_log(:).';          % 1 x p
    names  = GN.forMod_param_names;                      % 1 x p cell
    contrib_var = (sens.^2) .* var_b;                    % n_vars x p (log-space variances)

    for s = 1:n_src
        cols = strcmp(names, source_fields{s});
        var_s = sum(contrib_var(:, cols), 2);            % combined variance from this source
        frac_by_source(:, s, k) = sqrt(var_s);
        abs_by_source(:, s, k)  = x_hat .* sqrt(var_s);
    end

    % ---- sanity check: the 3 sources should sum (in variance) to diag(S_f) ----
    summed_var   = sum(frac_by_source(:, :, k).^2, 2);   % n_vars x 1
    decomp_residual(:, k) = summed_var - diag(GN.forward_model_parameter_cov);

    n_used = n_used + 1;
end

fprintf('Used %d pixels; skipped %d (empty/masked) and %d (missing decomposition fields).\n', ...
    n_used, n_skipped_empty, n_skipped_fields);

if n_used == 0
    error([newline, 'No usable pixels. If many were skipped for "missing fields", the ', ...
        'retrievals were run before the decomposition outputs were added to ', ...
        'calc_retrieval_gauss_newton_EMIT_ver4_log_forMo_uncert_perPixel.m.', newline])
end

% report the worst decomposition residual as a fraction of the total (should be ~0)
max_rel_resid = max(abs(decomp_residual(:)) ./ max(frac_fm_total(:).^2, eps), [], 'omitnan');
fprintf('Max relative decomposition residual (should be ~1e-10): %.2e\n', max_rel_resid);


%% ---- aggregate across pixels: median, 16th/84th percentiles, mean ----
% helper anonymous functions over the pixel dimension
med2  = @(A) median(A, 2, 'omitnan');
mean2 = @(A) mean(A, 2, 'omitnan');
p2    = @(A) prctile(A, [16 84], 2);                % n_vars x 2
med3  = @(A) median(A, 3, 'omitnan');
mean3 = @(A) mean(A, 3, 'omitnan');
p3    = @(A) prctile(A, [16 84], 3);                % n_vars x n_src x 2

budget = struct();
budget.n_pixels       = n_used;
budget.var_fields     = var_fields;
budget.source_fields  = source_fields;
budget.central_stat   = central_stat;

% --- fractional (relative) ---
budget.frac.total        = struct('median', med2(frac_total),    'mean', mean2(frac_total),    'p16_84', p2(frac_total));
budget.frac.fm_total     = struct('median', med2(frac_fm_total), 'mean', mean2(frac_fm_total), 'p16_84', p2(frac_fm_total));
budget.frac.by_source    = struct('median', med3(frac_by_source),'mean', mean3(frac_by_source),'p16_84', p3(frac_by_source));

% --- absolute (linear retrieval units) ---
budget.abs.total         = struct('median', med2(abs_total),     'mean', mean2(abs_total),     'p16_84', p2(abs_total));
budget.abs.fm_total      = struct('median', med2(abs_fm_total),  'mean', mean2(abs_fm_total),  'p16_84', p2(abs_fm_total));
budget.abs.by_source     = struct('median', med3(abs_by_source), 'mean', mean3(abs_by_source), 'p16_84', p3(abs_by_source));

budget.xhat              = struct('median', med2(xhat_all),      'mean', mean2(xhat_all),      'p16_84', p2(xhat_all));


%% ---- print the error budget ----
% pick central statistic for printing
if strcmp(central_stat, 'median')
    ctr_frac_total = budget.frac.total.median;     ctr_frac_src = budget.frac.by_source.median;
    ctr_frac_fm    = budget.frac.fm_total.median;
    ctr_abs_total  = budget.abs.total.median;      ctr_abs_src  = budget.abs.by_source.median;
    ctr_abs_fm     = budget.abs.fm_total.median;
else
    ctr_frac_total = budget.frac.total.mean;       ctr_frac_src = budget.frac.by_source.mean;
    ctr_frac_fm    = budget.frac.fm_total.mean;
    ctr_abs_total  = budget.abs.total.mean;        ctr_abs_src  = budget.abs.by_source.mean;
    ctr_abs_fm     = budget.abs.fm_total.mean;
end
plo_frac = budget.frac.by_source.p16_84(:, :, 1);  phi_frac = budget.frac.by_source.p16_84(:, :, 2);

fprintf('\n=====================================================================\n');
fprintf(' EMIT-Aqua retrieval error budget  (%s over %d pixels; [16th, 84th] pct)\n', central_stat, n_used);
fprintf('=====================================================================\n');
for i = 1:n_vars
    fprintf('\n%-6s   (retrieved %s = %.3g)\n', var_fields{i}, var_labels{i}, budget.xhat.median(i));
    fprintf('   %-26s : %6.2f%%   (abs %.3g)\n', 'TOTAL posterior',        100*ctr_frac_total(i), ctr_abs_total(i));
    fprintf('   %-26s : %6.2f%%   (abs %.3g)\n', 'forward-model (all)',    100*ctr_frac_fm(i),    ctr_abs_fm(i));
    for s = 1:n_src
        fprintf('     - %-22s : %6.2f%%  [%5.2f, %5.2f]   (abs %.3g)\n', ...
            source_labels{s}, 100*ctr_frac_src(i, s), 100*plo_frac(i, s), 100*phi_frac(i, s), ctr_abs_src(i, s));
    end
end
fprintf('\n(Per-source values are FRACTIONAL 1-sigma; "abs" is in the retrieved variable''s units.)\n');
fprintf('(Sanity: the 3 source variances sum to the forward-model total variance.)\n\n');


%% ---- optional figure: grouped bar chart of the fractional error budget ----
if plot_figures == true

    figure;
    bar_data = 100 * ctr_frac_src;                 % n_vars x n_src, percent
    b = bar(bar_data, 'grouped');                  % one Bar object per source (1 x n_src)
    for s = 1:n_src
        b(s).DisplayName = source_labels{s};
    end
    hold on

    % overlay the total forward-model and total posterior fractional uncertainty
    x_centers = (1:n_vars);
    h_fm  = plot(x_centers, 100*ctr_frac_fm,    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
        'LineStyle', 'none', 'DisplayName', 'forward-model (all)');
    h_tot = plot(x_centers, 100*ctr_frac_total, 'rd', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
        'LineStyle', 'none', 'DisplayName', 'total posterior');

    set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_labels, 'TickLabelInterpreter', 'tex');
    ylabel('Fractional 1\sigma uncertainty (%)', 'Interpreter', 'tex');
    title(sprintf('EMIT-Aqua forward-model error budget (%s of %d pixels)', central_stat, n_used));
    legend([b, h_fm, h_tot], 'Location', 'best', 'Interpreter', 'tex');
    grid on; grid minor
    set(gcf, 'Position', [100 100 850 500]);
end


%% ---- save aggregated results ----
if save_results == true
    out_name = fullfile(save_directory, ...
        ['EMIT_Aqua_forMod_error_budget_', char(datetime('today')), '.mat']);
    save(out_name, 'budget', 'retrieval_directory', 'n_used', ...
        'frac_total', 'frac_fm_total', 'frac_by_source', ...
        'abs_total', 'abs_fm_total', 'abs_by_source', 'xhat_all');
    fprintf('Saved error-budget results to:\n  %s\n', out_name);
end
