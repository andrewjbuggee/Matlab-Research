%% Backfill the forward-model-parameter error decomposition into the HySICS synthetic retrievals
%
% WHY THIS SCRIPT EXISTS
% ----------------------
% The paper-2 synthetic HySICS retrievals in paper2_variableSweep/..._uncert_3
% were produced by calc_retrieval_gauss_newton_HySICS_ver4_log_forMo_uncert_2.m,
% which — unlike its EMIT sibling — stops at the DOF/information-content outputs
% and NEVER ran the forward-model error decomposition block (EMIT lines 1214-1294).
% As a result each saved GN_outputs contains posterior_cov_log, Jacobian_final_log
% and retrieval, but NOT the pieces the error budget needs:
%       forward_model_parameter_cov, forMod_sensitivity_log,
%       forMod_param_variance_log,   forMod_param_names.
%
% The missing ingredient is the forward-model Jacobian K_b = d ln(R)/d ln(b),
% evaluated at the converged state. K_b CANNOT be reconstructed algebraically from
% the saved covariances (posterior_cov_log^-1 - S_a^-1 only recovers the 4x4
% product J'S_e^-1 J, which is not enough to invert back to the 636x636 S_e or to
% the 636x21 K_b). So this script RE-COMPUTES K_b once, at each saved converged
% retrieval state, by calling the identical forward-model Jacobian routine the
% retrieval used, then forms the Rodgers (2000) Eq. 3.18 decomposition and APPENDS
% the four fields back into each .mat file so the analysis script can read them.
%
% This re-runs libRadtran (the K_b routine perturbs all 21 forward-model
% parameters), so it must run where libRadtran AND the simulated-spectra files
% live — i.e. Alpine (curc). It is far cheaper than re-running the full
% Gauss-Newton retrievals (one K_b evaluation per file, not one per iteration).
%
% WHAT IT COMPUTES  (mirrors calc_retrieval_..._EMIT_..._perPixel.m lines 1217-1251)
%   G_log        = posterior_cov_log * K' * S_e^-1                 (gain, log space)
%   S_f          = G_log * K_b * S_b * K_b' * G_log'               (Eq. 3.18)
%   sensitivity  = G_log * K_b                                     (d ln x_hat / d ln b)
%   param_var    = diag(S_b)                                       (Var(ln b_j))
%   param_names  = [ re_profile x n_layers , cloud_top_height ]
% with K = Jacobian_final_log (saved) and S_e = S_y + K_b*S_b*K_b'.
%
% IDEMPOTENT: files that already carry the decomposition fields are skipped
% unless force_recompute = true.
%
% By Andrew John Buggee
%%

clear variables

% =============================== configuration ===============================
% Directory holding the synthetic retrievals to backfill. Set per machine; the
% curc branch is the one that matters because that is where libRadtran lives.
which_computer = whatComputer();

switch which_computer
    case 'curc'
        retrieval_directory = ['/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/', ...
            'Droplet_profile_retrievals/paper2_variableSweep/', ...
            'full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
    case 'andrewbuggee'
        retrieval_directory = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/', ...
            'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
    case 'anbu8374'
        retrieval_directory = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/', ...
            'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];
    otherwise
        error([newline, 'Define retrieval_directory for computer "', which_computer, '".', newline])
end

% Scratch/INP folder-extension number handed to define_folderPaths_for_HySICS.
% Pick one that will not collide with a running retrieval on the same machine.
folder_extension_number = 999;

force_recompute = false;   % true -> re-backfill even files that already have the fields
save_backfill   = true;    % true -> append the four fields into each .mat via -append
% =============================================================================


%% ---- build machine-valid folder paths and a parallel pool ----
% The stale folder_paths saved in each file point at scratch dirs from the
% original SLURM run; rebuild fresh, valid, writable paths for THIS session.
folder_paths = define_folderPaths_for_HySICS(folder_extension_number);
start_parallel_pool(folder_paths.which_computer);

% jacobian_barPlot_flag is a plotting switch inside the K_b routine; keep it off.
jacobian_barPlot_flag = false;


%% ---- load the instrument spectral response once (fixed across all pixels) ----
% spec_response is an instrument property, identical for every 636-band HySICS
% file, and is NOT stored in the saved retrievals — grab it from any simulated
% measurement file in the simulated-spectra folder.
spec_files = dir(fullfile(folder_paths.HySICS_simulated_spectra, '*.mat'));
if isempty(spec_files)
    error([newline, 'No simulated-spectra .mat found in:', newline, ...
        folder_paths.HySICS_simulated_spectra, newline, ...
        'spec_response cannot be loaded — point HySICS_simulated_spectra at the sim folder.', newline])
end
sm = load(fullfile(spec_files(1).folder, spec_files(1).name), 'spec_response');
if ~isfield(sm, 'spec_response') || ~isfield(sm.spec_response, 'value')
    error([newline, 'spec_response.value not found in ', spec_files(1).name, newline])
end
spec_response = sm.spec_response.value;
fprintf('Loaded spec_response (%d bands) from %s\n', size(spec_response, 1), spec_files(1).name);


%% ---- find the retrieval files ----
file_list = dir(fullfile(retrieval_directory, 'dropletRetrieval_HySICS_*.mat'));
if isempty(file_list)
    error([newline, 'No dropletRetrieval_HySICS_*.mat files in:', newline, retrieval_directory, newline])
end
n_files = numel(file_list);
fprintf('Found %d retrieval files to backfill in\n  %s\n\n', n_files, retrieval_directory);

decomp_fields = {'forward_model_parameter_cov', 'forMod_sensitivity_log', ...
    'forMod_param_variance_log', 'forMod_param_names'};

n_done = 0; n_skipped = 0; n_failed = 0;
max_rel_post_resid = 0;   % running worst mismatch between recomputed & saved posterior


%% ---- loop over pixels ----
for k = 1:n_files

    fpath = fullfile(file_list(k).folder, file_list(k).name);
    fprintf('[%d/%d] %s\n', k, n_files, file_list(k).name);

    try
        S = load(fpath, 'GN_outputs', 'GN_inputs');
        if ~isfield(S, 'GN_outputs') || ~isfield(S, 'GN_inputs')
            warning('  missing GN_outputs/GN_inputs — skipped.'); n_failed = n_failed + 1; continue
        end
        GN_outputs = S.GN_outputs;
        GN_inputs  = S.GN_inputs;

        % skip already-backfilled files unless forced
        if all(isfield(GN_outputs, decomp_fields)) && ~force_recompute
            fprintf('  already has decomposition fields — skipped.\n'); n_skipped = n_skipped + 1; continue
        end

        % ---- inputs to the K_b routine, exactly as the retrieval called it ----
        state_vector           = GN_outputs.retrieval(:, end);          % [r_top; r_bot; tau_c; acpw]
        measurement_estimate_ln = GN_outputs.new_measurement_estimate;  % ln-space meas estimate at convergence
        GN_inputs.which_computer = folder_paths.which_computer;         % override stale machine tag
        GN_inputs.spec_response  = spec_response;                       % some internal branches read this

        % ---- re-compute the forward-model Jacobian K_b (636 x 21) ----
        % Columns: [ r_e profile per layer (cloud-top -> base) , cloud top height ].
        % GN_inputs.model.forward_model.re.mean{end} was appended to the converged
        % profile at retrieval time and is stored in the saved GN_inputs, so the
        % routine reads the correct converged profile.
        K_b = compute_forMod_jacobian_HySICS_log_reProf_cloudTopHeight(state_vector, ...
            measurement_estimate_ln, GN_inputs, spec_response, jacobian_barPlot_flag, folder_paths);

        if any(~isfinite(K_b), 'all')
            warning('  K_b contains non-finite entries — skipped.'); n_failed = n_failed + 1; continue
        end

        % ---- covariance ingredients (all in log space) ----
        S_b   = GN_inputs.model.forward_model.covariance;   % 21x21 diagonal, Var(ln b)
        S_y   = GN_inputs.measurement.covariance;           % 636x636 measurement cov
        S_a   = GN_inputs.model.covariance;                 % 4x4 prior (a priori) cov
        J_log = GN_outputs.Jacobian_final_log;              % 636x4, saved converged Jacobian

        % effective measurement covariance with the forward-model term folded in
        %   S_e = S_y + K_b * S_b * K_b'
        S_e = S_y + (K_b * S_b * K_b');

        % re-derive the posterior covariance from the recomputed K_b and compare to
        % the saved one — a consistency check that K_b was reproduced correctly.
        Se_inv = S_e^(-1);
        posterior_cov_log_recomp = ((J_log' * Se_inv * J_log) + S_a^(-1))^(-1);
        rel_post_resid = max(abs(posterior_cov_log_recomp(:) - GN_outputs.posterior_cov_log(:))) ...
            / max(abs(GN_outputs.posterior_cov_log(:)));
        max_rel_post_resid = max(max_rel_post_resid, rel_post_resid);
        if rel_post_resid > 1e-3
            warning('  recomputed posterior differs from saved by %.2e (K_b reproduction imperfect).', rel_post_resid);
        end

        % ---- Rodgers (2000) Eq. 3.18 forward-model error covariance (log space) ----
        % use the self-consistent recomputed posterior so G and S_f share one K_b.
        G_log = posterior_cov_log_recomp * J_log' * Se_inv;                 % gain, 4x636
        forward_model_parameter_cov = G_log * K_b * S_b * K_b' * G_log';    % 4x4, S_f

        % ---- per-parameter decomposition ingredients ----
        forMod_sensitivity_log   = G_log * K_b;             % 4x21, d ln x_hat / d ln b_j
        forMod_param_variance_log = diag(S_b);              % 21x1, Var(ln b_j)

        % physical grouping of the 21 columns: first n_layers are the r_e profile,
        % the last one is cloud top height (matches the K_b column assembly).
        n_layers = GN_inputs.RT.n_layers;
        if numel(forMod_param_variance_log) ~= n_layers + 1
            warning('  expected %d forward-model params (n_layers+1) but S_b has %d.', ...
                n_layers + 1, numel(forMod_param_variance_log));
        end
        forMod_param_names = [repmat({'re_profile'}, 1, numel(forMod_param_variance_log) - 1), ...
            {'cloud_top_height'}];

        % ---- write the fields back ----
        GN_outputs.forward_model_parameter_cov = forward_model_parameter_cov;
        GN_outputs.forMod_sensitivity_log      = forMod_sensitivity_log;
        GN_outputs.forMod_param_variance_log   = forMod_param_variance_log;
        GN_outputs.forMod_param_names          = forMod_param_names;

        if save_backfill
            save(fpath, 'GN_outputs', '-append');
        end

        fm_frac = sqrt(diag(forward_model_parameter_cov));   % 1-sigma fractional FM uncertainty
        fprintf('  OK. FM 1sigma [r_top r_bot tau_c acpw] = [%.3f %.3f %.3f %.3f]  (posterior match %.1e)\n', ...
            fm_frac(1), fm_frac(2), fm_frac(3), fm_frac(4), rel_post_resid);
        n_done = n_done + 1;

    catch ME
        warning('backfill:pixelFailed', '  FAILED: %s', ME.message);
        n_failed = n_failed + 1;
    end
end

fprintf('\n==================== backfill complete ====================\n');
fprintf('  backfilled : %d\n  skipped    : %d (already had fields)\n  failed     : %d\n', ...
    n_done, n_skipped, n_failed);
fprintf('  worst recomputed-vs-saved posterior mismatch: %.2e (should be << 1e-3)\n\n', max_rel_post_resid);
