function summary = run_synthetic_spectra_chunk(input_file, cloud_ids, folder_paths, output_dir, settings)
%% HySICS spectra for a chunk of synthetic clouds, parallel over clouds, one output file
%
% Runs hysics_refl_from_synthetic_NN_inputs_ver2 for every cloud in
% cloud_ids with parfor over CLOUDS (each worker computes whole clouds), and
% writes ONE .mat file for the chunk instead of one file per cloud.
%
% Why
%   * Parallel over clouds: each cloud is a handful of independent uvspec runs,
%     so whole clouds are the natural unit of work for a worker.
%   * One file per chunk: per-cloud files were ~2.6 MB each, almost all of it
%     a copy of the 636 x 401 response-function table; 100,000 of them is
%     ~260 GB in 100,000 files. CURC's filesystem guidance: "If you can choose
%     between one file and 100 files for storing the same amount of data, you
%     will see far better performance from the single large file"
%     (curc.readthedocs.io, compute/filesystems.html).
%   * Checkpointing: clouds are processed in blocks and the file is rewritten
%     after every block, so a task that hits its wall-time limit loses at most
%     one block. On restart, clouds already marked done in an existing file are
%     skipped.
%
% INPUTS
%   input_file    training-input NetCDF (07_build_training_inputs.py)
%   cloud_ids     vector of 1-based cloud indices
%   folder_paths  from define_folderPaths_for_HySICS
%   output_dir    folder for the chunk file (trailing '/')
%   settings      (optional) struct:
%                   .n_segments   uvspec runs per cloud (default 4)
%                   .fwhm_nm      HySICS FWHM (default 6)
%                   .block_size   clouds per checkpoint (default 2 x pool size)
%                   .file_tag     text added to the file name (default '')
%                   .verbose      per-segment printing (default false)
%
% OUTPUT FILE  <output_dir>synthetic_spectra_HySICS_636bands_clouds<first>-<last><tag>.mat
%   cloud_idx (n x 1), done (n x 1 logical), err_msg (n x 1 cell),
%   Refl_model (n x 636) clean reflectance (1/sr; NaN where not done),
%   sza saz vza vaz (n x 1, deg), re_um lwc_g_per_m3 z_km (n x 7, top -> base),
%   tau_c LWP_g_per_m2 alpha (n x 1), T_K vapor_cm3 (n x 37), pressure_hPa (37 x 1),
%   mie_table (n x 1 cell), t_cloud_s (n x 1), band_center_wavelength_nm,
%   band_fwhm_nm, band_window_nm (636 x 2), input_file, settings, created
%
% By Andrew John Buggee; written with Claude (assistant), 14 Nov 2026
%%

if nargin < 5 || isempty(settings), settings = struct(); end
if ~isfield(settings, 'n_segments'), settings.n_segments = 4;  end
if ~isfield(settings, 'fwhm_nm'),    settings.fwhm_nm    = 6;  end
if ~isfield(settings, 'file_tag'),   settings.file_tag   = ''; end
if ~isfield(settings, 'verbose'),    settings.verbose    = false; end
settings.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';

which_computer = whatComputer();
cloud_ids = cloud_ids(:);
n = numel(cloud_ids);

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
out_file = sprintf('%ssynthetic_spectra_HySICS_636bands_clouds%d-%d%s.mat', output_dir, ...
    cloud_ids(1), cloud_ids(end), settings.file_tag);

% ---- once per chunk: response functions and solar flux ----------------------
spec_response = create_HySICS_specResponse_ver2((1:636)', settings.source_file, ...
    which_computer, settings.fwhm_nm);
[settings.source_flux, settings.source_wavelength] = read_solar_flux_file( ...
    [min(spec_response.window, [], 'all'), max(spec_response.window, [], 'all')], ...
    settings.source_file);

% ---- results arrays (resume from an existing chunk file if present) ----------
R = init_results(n, cloud_ids);
if isfile(out_file)
    prev = load(out_file);
    if isfield(prev, 'cloud_idx') && isequal(prev.cloud_idx(:), cloud_ids)
        R = copy_fields(R, prev);
        fprintf('Resuming %s: %d of %d clouds already done\n', out_file, nnz(R.done), n);
    end
end

% ---- parallel pool (started by the batch script; start one if absent) -------
pool = gcp('nocreate');
if isempty(pool)
    start_parallel_pool(which_computer);
    pool = gcp('nocreate');
end
if ~isfield(settings, 'block_size') || isempty(settings.block_size)
    settings.block_size = 2 * pool.NumWorkers;
end

todo = find(~R.done);
fprintf('Chunk %d-%d: %d clouds to run, %d workers, blocks of %d, %d segments per cloud\n', ...
    cloud_ids(1), cloud_ids(end), numel(todo), pool.NumWorkers, settings.block_size, ...
    settings.n_segments);

drv_settings = rmfield(settings, {'fwhm_nm', 'file_tag', 'block_size'});

for b0 = 1:settings.block_size:numel(todo)

    blk  = todo(b0:min(b0 + settings.block_size - 1, numel(todo)));
    nb   = numel(blk);
    outs = cell(nb, 1);
    errs = cell(nb, 1);
    t_blk = tic;

    parfor j = 1:nb
        try
            outs{j} = hysics_refl_from_synthetic_NN_inputs_ver2(input_file, cloud_ids(blk(j)), ...
                folder_paths, spec_response, drv_settings);
            errs{j} = '';
        catch ME
            errs{j} = sprintf('%s: %s', ME.identifier, ME.message);
        end
    end

    for j = 1:nb
        k = blk(j);
        if isempty(errs{j})
            R = store_cloud(R, k, outs{j});
        else
            R.err_msg{k} = errs{j};
            fprintf('[CLOUD %d FAILED] %s\n', cloud_ids(k), errs{j});
        end
    end

    save_chunk(out_file, R, spec_response, input_file, settings);
    fprintf('[%s] block of %d clouds in %.1f s; %d/%d done\n', char(datetime('now')), nb, ...
        toc(t_blk), nnz(R.done), n);
end

summary.out_file = out_file;
summary.n_done   = nnz(R.done);
summary.n_failed = nnz(~R.done);
summary.median_t_cloud_s = median(R.t_cloud_s(R.done));
fprintf('\n=== CHUNK SUMMARY: %d done, %d failed, median %.1f s per cloud -> %s ===\n', ...
    summary.n_done, summary.n_failed, summary.median_t_cloud_s, out_file);

end


% =============================================================================
function R = init_results(n, cloud_ids)
R.cloud_idx    = cloud_ids;
R.done         = false(n, 1);
R.err_msg      = repmat({''}, n, 1);
R.Refl_model   = nan(n, 636);
R.sza = nan(n, 1); R.saz = nan(n, 1); R.vza = nan(n, 1); R.vaz = nan(n, 1);
R.re_um        = nan(n, 7);
R.lwc_g_per_m3 = nan(n, 7);
R.z_km         = nan(n, 7);
R.tau_c = nan(n, 1); R.LWP_g_per_m2 = nan(n, 1); R.alpha = nan(n, 1);
R.T_K          = nan(n, 37);
R.vapor_cm3    = nan(n, 37);
R.pressure_hPa = nan(37, 1);
R.mie_table    = repmat({''}, n, 1);
R.t_cloud_s    = nan(n, 1);
end


function R = copy_fields(R, prev)
% copy every results field present in a previous chunk file
f = fieldnames(R);
for i = 1:numel(f)
    if isfield(prev, f{i}), R.(f{i}) = prev.(f{i}); end
end
end


function R = store_cloud(R, k, o)
R.done(k)            = true;
R.err_msg{k}         = '';
R.Refl_model(k, :)   = o.refl(:)';
R.sza(k) = o.sza; R.saz(k) = o.saz; R.vza(k) = o.vza; R.vaz(k) = o.vaz;
R.re_um(k, :)        = o.re_um(:)';
R.lwc_g_per_m3(k, :) = o.lwc_g_per_m3(:)';
R.z_km(k, :)         = o.z_km(:)';
R.tau_c(k) = o.tau_c; R.LWP_g_per_m2(k) = o.LWP_g_per_m2; R.alpha(k) = o.alpha;
R.T_K(k, :)          = o.T_K(:)';
R.vapor_cm3(k, :)    = o.vapor_cm3(:)';
R.pressure_hPa       = o.pressure_hPa(:);
R.mie_table{k}       = o.mie_table;
R.t_cloud_s(k)       = o.t_total_s;
end


function save_chunk(out_file, R, spec_response, input_file, settings)
% write to a temporary name and rename, so an interrupted save never leaves a
% truncated chunk file behind
S = R;
S.band_center_wavelength_nm = spec_response.center_wavelength;
S.band_fwhm_nm              = spec_response.fwhm;
S.band_window_nm            = spec_response.window;
S.input_file                = input_file;
S.settings                  = rmfield(settings, {'source_flux', 'source_wavelength'});
S.created                   = char(datetime('now'));
tmp = [out_file(1:end-4), '_partial.mat'];
save(tmp, '-struct', 'S', '-v7');
movefile(tmp, out_file, 'f');
end
