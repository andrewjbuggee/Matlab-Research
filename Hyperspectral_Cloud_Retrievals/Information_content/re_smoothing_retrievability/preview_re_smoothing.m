%% Preview the vertical smoothing of in-situ effective-radius profiles
%
% Run this LOCALLY (on your Mac) BEFORE submitting the supercomputer jobs.
% It loads the VOCALS-REx and ORACLES in-situ ensemble profiles, applies the
% moving-average windows defined in define_re_smoothing_windows.m, and lets
% you see how aggressively each window removes the high-frequency vertical
% structure in r_e(z). Tweak the windows in define_re_smoothing_windows.m and
% re-run until you are happy, then submit the batch jobs.
%
% It produces:
%   (1) Per-campaign example panels: raw r_e(z) vs each smoothed version, for
%       a few profiles spanning low / medium / high vertical variability.
%   (2) An aggregate summary: across ALL profiles, the fractional reduction in
%       the vertical standard deviation of r_e for each window. This is a
%       single-number measure of "how much fluctuation each window removes."
%
% Nothing here runs libRadtran. It is purely a visualization of the smoothing.
%
% By Andrew John Buggee

clear variables; close all

%% ---------------------------- CONFIG --------------------------------------

% smoothing windows (meters) - edit define_re_smoothing_windows.m to change
windows_m = define_re_smoothing_windows();

% how many example profiles to show per campaign
n_examples_per_campaign = 3;

% save the figures?
save_figures = false;

%% ----------------------- locate the data ----------------------------------

which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374')
    base = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/';
elseif strcmp(which_computer, 'andrewbuggee')
    base = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/';
else
    error('This preview is meant to be run on your local machine.')
end

vocals_file = [base, 'VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/', ...
    'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_', ...
    'drizzleLWP-threshold_5_10-Nov-2025.mat'];

oracles_file = [base, 'ORACLES/oracles_data/', ...
    'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_', ...
    'no_rEff_greaterThan50_microns_16-Mar-2026.mat'];

fig_folder = [base, 'Information_content/re_smoothing_retrievability/preview_figures/'];
if save_figures && ~exist(fig_folder, 'dir'); mkdir(fig_folder); end

campaigns = { 'VOCALS-REx', vocals_file; ...
              'ORACLES',    oracles_file };

%% ----------------- loop over campaigns: examples + stats -------------------

ax_lbl = 18; 
title_lbl = 18;

for cc = 1:size(campaigns, 1)

    campaign_name = campaigns{cc, 1};
    ds = load(campaigns{cc, 2});
    ep = ds.ensemble_profiles;
    n_prof = numel(ep);

    % --- pull each profile's cleaned (re, z) and its vertical variability ---
    re_all = cell(n_prof, 1);
    z_all  = cell(n_prof, 1);
    std_raw = nan(n_prof, 1);                       % vertical std of raw re
    std_frac_remaining = nan(n_prof, numel(windows_m));   % std_smooth / std_raw

    for pp = 1:n_prof

        [re_um, z_m] = local_get_clean_profile(ep{pp});
        if numel(re_um) < 3; continue; end

        re_all{pp} = re_um;
        z_all{pp}  = z_m;
        std_raw(pp) = std(re_um);

        for ww = 1:numel(windows_m)
            re_s = smooth_re_profile(re_um, z_m, windows_m(ww));
            if std_raw(pp) > 0
                std_frac_remaining(pp, ww) = std(re_s) / std_raw(pp);
            end
        end
    end

    % --- choose example profiles spanning low/med/high variability ---
    valid = find(~isnan(std_raw));
    [~, order] = sort(std_raw(valid), 'ascend');
    valid_sorted = valid(order);
    if numel(valid_sorted) >= n_examples_per_campaign
        pick_pos = round(linspace(1, numel(valid_sorted), n_examples_per_campaign));
        examples = valid_sorted(pick_pos);
    else
        examples = valid_sorted;
    end

    % ------------------------- plot the examples -------------------------
    figure('Name', [campaign_name, ' - r_e smoothing examples'], ...
        'Position', [80, 80, 380*numel(examples), 560]);
    tl = tiledlayout(1, numel(examples), 'TileSpacing', 'compact', 'Padding', 'compact');

    cmap = parula(numel(windows_m) + 1);

    for ee = 1:numel(examples)
        pp = examples(ee);
        re_um = re_all{pp};
        z_m   = z_all{pp};

        nexttile; hold on
        % raw
        plot(re_um, z_m, '-o', 'Color', 'k', 'LineWidth', 2, ...
            'MarkerSize', 3, 'DisplayName', 'raw');
        % smoothed versions
        for ww = 1:numel(windows_m)
            re_s = smooth_re_profile(re_um, z_m, windows_m(ww));
            plot(re_s, z_m, '-', 'Color', cmap(ww+1, :), 'LineWidth', 2, ...
                'DisplayName', sprintf('%d m', windows_m(ww)));
        end
        grid on; box on
        xlabel('$r_e$  $[\mu m]$', 'FontSize', ax_lbl, 'Interpreter','latex')
        if ee == 1; ylabel('Altitude  $[m]$', 'FontSize', ax_lbl, 'Interpreter','latex'); end
        title(sprintf('profile %d  ($\\sigma_{r_e}=%.2f \\mu m$)', pp, std(re_um)),...
            'FontSize', ax_lbl, 'Interpreter','latex')
        if ee == numel(examples); legend('Location', 'best'); end
    end
    title(tl, [campaign_name, ': raw vs vertically smoothed $r_e(z)$'], ...
        'FontWeight', 'bold', 'FontSize', ax_lbl, 'Interpreter','latex')

    if save_figures
        saveas(gcf, [fig_folder, 'examples_', campaign_name, '.png']);
    end

    % --------------- aggregate aggressiveness across profiles ---------------
    % fractional REDUCTION in vertical std of r_e: 0 = no change, 1 = flat
    frac_reduction = 1 - std_frac_remaining;   % n_prof x numel(windows)

    figure('Name', [campaign_name, ' - smoothing aggressiveness'], ...
        'Position', [120, 120, 640, 460]);
    boxplot(frac_reduction, 'Labels', compose('%d m', windows_m));
    grid on; box on
    ylim([0, 1])
    xlabel('moving-average window', 'FontSize', ax_lbl, 'Interpreter','latex')
    ylabel('fractional reduction in vertical $\sigma_{r_e}$', 'FontSize', ax_lbl, 'Interpreter','latex')
    title({[campaign_name, ': how much vertical $r_e$ variability each window removes'], ...
        sprintf('(across %d profiles; 1 = profile flattened to its mean)', numel(valid))},...
        'FontSize', ax_lbl, 'Interpreter','latex')

    if save_figures
        saveas(gcf, [fig_folder, 'aggressiveness_', campaign_name, '.png']);
    end

    % ------------------------- print a quick summary -------------------------
    fprintf('\n==== %s (%d valid profiles) ====\n', campaign_name, numel(valid));
    fprintf('  window [m] | median frac. std removed | %% profiles >90%% flattened\n');
    for ww = 1:numel(windows_m)
        med_removed = median(frac_reduction(valid, ww));
        pct_flat = 100 * mean(frac_reduction(valid, ww) > 0.90);
        fprintf('  %8d   |        %5.1f%%            |        %5.1f%%\n', ...
            windows_m(ww), 100*med_removed, pct_flat);
    end
end

fprintf('\nDone. If the windows look right, submit the batch jobs.\n');
fprintf('If not, edit define_re_smoothing_windows.m and re-run this script.\n');


%% ===================== local helper function ==============================

function [re_um, z_m] = local_get_clean_profile(profile)
% Extract and clean a single in-situ (re, z) profile the SAME way the
% reflectance calculation does: pick the effective-radius field, drop levels
% with r_e >= 50 um (outside the custom Mie table) or r_e < 1 um (cloud-edge
% artifacts / below the Mie table), and return altitude in meters.

    % --- effective radius field (microns) ---
    if isfield(profile, 're')
        re_um = profile.re(:);
    elseif isfield(profile, 're_CDP')
        re_um = profile.re_CDP(:);
    else
        re_um = [];
        z_m   = [];
        return
    end

    % --- altitude (meters) ---
    z_m = profile.altitude(:);

    % --- drop NaNs ---
    bad = isnan(re_um) | isnan(z_m);
    re_um(bad) = []; z_m(bad) = [];

    % --- drop out-of-range droplet sizes (mirror the RT-side cleaning) ---
    drop = re_um >= 50 | re_um < 1;
    re_um(drop) = []; z_m(drop) = [];
end
