#!/bin/bash
#
# PRODUCTION Slurm job array for the 12-Nov-2026 synthetic-cloud RT training
# set (100,000 clouds), ver2 pipeline (14 Nov 2026):
#   run_synthetic_spectra_chunk.m -> hysics_refl_from_synthetic_NN_inputs_ver2.m
# Each cloud: 4 uvspec runs over contiguous groups of bands instead of 636
# per-band runs (identical reflectances, ~11x less CPU); parfor over clouds;
# one output .mat per task, rewritten after every block of clouds so a task
# that times out keeps what it finished, and resubmitting the same task
# resumes from that file.
#
# Cloud index range covered by task t (1-indexed):
#   start = CLOUD_OFFSET + (t-1)*CHUNK_SIZE + 1
#   end   = min(CLOUD_OFFSET + t*CHUNK_SIZE, N_TOTAL)
# With CHUNK_SIZE = 600: ceil(100000/600) = 167 tasks, one array (CURC limit
# 1,000 tasks), --array=1-167%40.
#
# CHUNK SIZE vs WALL TIME: set CHUNK_SIZE from the test's median time per
# cloud on Alpine (t_cloud, printed as "median N s per cloud", measured on
# 12 workers in the test):  CHUNK_SIZE ~ 0.5 x CPUS x WALL_s / t_cloud  (factor 0.5 =
# headroom). The default 600 assumes t_cloud ~ 700 s on 32 CPUs in 8 h:
# 600 x 700 / 32 = 3.6 h expected. Re-derive the --array upper bound after
# changing it: ceil(N_TOTAL / CHUNK_SIZE).
#
# Memory: the largest segment run (332-856 nm) peaked at 3.1 GB on the laptop
# (a single full-range run: 5.0 GB), plus ~1 GB per MATLAB worker, so budget
# ~4 GB per worker: 32 x 4 + ~3 GB client = ~131 GB -> --mem=135G. That is
# above 32 x ~3.75 GB, so CURC bills ~36 cores (job-resources page: memory
# beyond the per-core share is charged as cores).
#
# BEFORE SUBMITTING (on Alpine):
#   1. module load slurm/alpine      (on login-ci4 this was needed for sinfo to
#                                     list acpu; the CURC batch-jobs page does not
#                                     mention it, so it may be harmless if absent)
#   2. git pull in /projects/anbu8374/Matlab-Research   (needs run_synthetic_spectra_chunk.m)
#   3. The input .nc is on /scratch (the test job already used it)
#   4. mkdir -p /projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch1
#      (Slurm cannot create the --chdir directory)
#
# OUTPUT RETENTION: /scratch/alpine files are "automatically removed 90 days
# after their initial creation" (CURC filesystems page, read 13 Nov 2026).
# Copy the finished training set off /scratch (to PetaLibrary or the laptop)
# well before then — the first clouds written are the first to expire.
#
# Partition / QoS: acpu accepts cpu-normal, cpu-long and cpu-testing (sbatch
# error message, 12 Nov 2026). cpu-normal with the partition's 1-day limit
# matches the old 24-h tasks; use cpu-long only if a chunk needs longer.
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --mem=135G              # ~4 GB per worker (3.1 GB uvspec + MATLAB); billed as ~36 cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=synthetic_spectra_ver2
#SBATCH --chdir=/projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch1
#SBATCH --output=create_meas_synthetic_NN_trainingData_12Nov_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_12Nov_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-167%40

# --- Chunk parameters (edit or override between batches) -----------------
CHUNK_SIZE=${CHUNK_SIZE_OVERRIDE:-600}        # clouds per array task
CLOUD_OFFSET=${CLOUD_OFFSET_OVERRIDE:-0}      # one array covers all clouds
N_SEGMENTS=${N_SEGMENTS_OVERRIDE:-4}          # uvspec runs per cloud
FILE_TAG=""
N_TOTAL=${N_TOTAL_OVERRIDE:-100000}           # length of the 'cloud' dimension in the input .nc
# -------------------------------------------------------------------------

# Modules
ml purge
ml gcc/11.2.0
ml netcdf/4.8.1
ml perl/5.36.0
ml texlive/2021

# libRadtran + GSL
export PATH=/projects/$USER/software/libRadtran-2.0.5/:$PATH
export PATH=/projects/$USER/software/libRadtran-2.0.5/data/:$PATH
export PATH=/projects/$USER/software/libRadtran-2.0.5/bin/:$PATH
export GSL_BIN=/projects/$USER/software/gsl-2.6/bin
export GSL_LIB=/projects/$USER/software/gsl-2.6/lib
export GSL_INC=/projects/$USER/software/gsl-2.6/include
export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH
export INSTALL_DIR=/projects/$USER/software/libRadtran-2.0.5
export PATH=$GSL_BIN:$PATH
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

cd /projects/anbu8374/ || exit 1

module load matlab/R2024b

# Input file (the .nc produced by 07_build_training_inputs.py on 12 Nov 2026)
input_file="/scratch/alpine/anbu8374/neural_network_training_data/training_inputs_jointMVN_N100000_L7_2026-11-12_ver2.nc"

# Output dir for the whole 12-Nov-2026 set (both batches write here).
# (Trailing slash required by the MATLAB driver.)
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/synthetic_dataSet_created_on_12_Nov_2026/"
mkdir -p "$output_dir"

# --- Pre-flight checks: fail in seconds, not after MATLAB starts -----------
code_dir="/projects/anbu8374/Matlab-Research"
echo "Matlab-Research commit: $(git -C "$code_dir" log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo 'unknown (not a git checkout?)')"
if [ -z "$(find "$code_dir" -name write_wc_file_mixed_layers.m -print -quit)" ]; then
    echo "ERROR: write_wc_file_mixed_layers.m not found under $code_dir. git pull first." >&2
    exit 1
fi
if [ ! -f "$input_file" ]; then
    echo "ERROR: input file not found: $input_file" >&2
    exit 1
fi
if [ -z "$(find "$code_dir" -name run_synthetic_spectra_chunk.m -print -quit)" ]; then
    echo "ERROR: run_synthetic_spectra_chunk.m not found under $code_dir. git pull first." >&2
    exit 1
fi
# Guard against the pre-azimuth-fix input file: the corrected
# training_inputs_jointMVN_N100000_L7_2026-11-12_ver2.nc (07, 12 Nov 2026,
# 22:04) has this MD5; the old wrong-azimuth _ver2 was d8e67d5b5db5157e9ff565fc958d0f1e.
EXPECTED_MD5=${EXPECTED_MD5_OVERRIDE:-f619a48cffea9631ea3899faa1166af4}
actual_md5=$(md5sum "$input_file" | cut -d' ' -f1)
if [ "$actual_md5" != "$EXPECTED_MD5" ]; then
    echo "ERROR: $input_file has MD5 $actual_md5, expected $EXPECTED_MD5 (azimuth-corrected file)." >&2
    exit 1
fi
echo "Input MD5 OK: $actual_md5"
# -------------------------------------------------------------------------

# Compute this task's cloud-index range, with clamping at N_TOTAL.
start_id=$(( CLOUD_OFFSET + (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
end_id=$((   CLOUD_OFFSET + SLURM_ARRAY_TASK_ID       * CHUNK_SIZE     ))
if [ "$end_id" -gt "$N_TOTAL" ]; then
    end_id=$N_TOTAL
fi
if [ "$start_id" -gt "$N_TOTAL" ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID}: start_id=$start_id is past N_TOTAL=$N_TOTAL. Nothing to do."
    exit 0
fi
n_clouds=$(( end_id - start_id + 1 ))
echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${n_clouds} clouds (ids ${start_id}..${end_id})"

# Per-task scratch on node-local /tmp instead of Lustre. This moves parpool's
# JobStorageLocation off the shared filesystem and avoids the multi-minute
# "Job Queued" parpool waits seen in earlier runs. Falls back to Lustre if
# /tmp is unwritable for any reason.
# CURC recommends $SLURM_SCRATCH: the job's own directory on the node's local
# SSD (~300 GB/node), deleted automatically when the job ends
# (curc.readthedocs.io, compute/filesystems.html, read 13 Nov 2026).
# Falls back to /tmp if it is not set.
TMPDIR_LOCAL="${SLURM_SCRATCH:-/tmp}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if mkdir -p "$TMPDIR_LOCAL" 2>/dev/null && [ -w "$TMPDIR_LOCAL" ]; then
    export TMPDIR="$TMPDIR_LOCAL"
    echo "TMPDIR (node-local, SLURM_SCRATCH=${SLURM_SCRATCH:-unset}): $TMPDIR"
else
    export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
    echo "TMPDIR (Lustre fallback): $TMPDIR"
fi

# Report headroom on the chosen TMPDIR so capacity problems show in the .out log.
echo "TMPDIR headroom:"
df -h "$TMPDIR"

# --- Stage the custom Mie tables on node-local disk (added 13 Nov 2026) ----
# Every uvspec run opens the cloud's custom Mie table (~130 MB each, 1-50 um
# family): 4 runs per cloud with the ver2 driver, 636 with ver1. CURC asks
# that /projects not carry intensive job I/O. Copy the 48 tables once per task to $SLURM_SCRATCH (local SSD, deleted
# at job end) and point MATLAB at them via CUSTOM_MIE_TABLES_DIR, which
# find_custom_mieTable_closest_to_alpha_profile.m honours. Source: a /scratch
# copy if one exists (Lustre is built for job I/O), otherwise /projects.
# If the copy is incomplete, fall back to the original tables with a warning.
MIE_PATTERN="wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_*.cdf"
MIE_N_EXPECTED=48
MIE_SCRATCH_SRC="/scratch/alpine/${USER}/custom_mieTables"
MIE_PROJECTS_SRC="/projects/${USER}/software/libRadtran-2.0.5/data/wc/custom_mieTables"
MIE_SRC=""
for cand in "$MIE_SCRATCH_SRC" "$MIE_PROJECTS_SRC"; do
    if [ "$(find "$cand" -maxdepth 1 -name "$MIE_PATTERN" 2>/dev/null | wc -l)" -ge "$MIE_N_EXPECTED" ]; then
        MIE_SRC="$cand"
        break
    fi
done
# Say loudly when the /scratch copy is gone (e.g. removed by the 90-day
# /scratch purge) so it can be re-copied; the job still runs from /projects.
if [ "$MIE_SRC" != "$MIE_SCRATCH_SRC" ]; then
    n_scratch=$(find "$MIE_SCRATCH_SRC" -maxdepth 1 -name "$MIE_PATTERN" 2>/dev/null | wc -l | tr -d ' ')
    msg="NOTICE: /scratch Mie-table copy is missing or incomplete (${n_scratch} of ${MIE_N_EXPECTED} files in ${MIE_SCRATCH_SRC}); staging from ${MIE_SRC:-nowhere} instead. Re-copy on a login node with:
    mkdir -p ${MIE_SCRATCH_SRC} && cp ${MIE_PROJECTS_SRC}/${MIE_PATTERN} ${MIE_SCRATCH_SRC}/"
    echo "$msg"          # .out log
    echo "$msg" >&2      # .err log, where it stands out
fi
MIE_LOCAL="${SLURM_SCRATCH:-$TMPDIR}/custom_mieTables"
if [ -n "$MIE_SRC" ] && mkdir -p "$MIE_LOCAL"; then
    t_copy=$(date +%s)
    echo "Mie tables: copying ${MIE_N_EXPECTED} files from $MIE_SRC to $MIE_LOCAL (started $(date '+%H:%M:%S')) ..."
    cp "$MIE_SRC"/$MIE_PATTERN "$MIE_LOCAL"/ 2>/dev/null
    n_copied=$(find "$MIE_LOCAL" -maxdepth 1 -name "$MIE_PATTERN" | wc -l)
    if [ "$n_copied" -ge "$MIE_N_EXPECTED" ]; then
        export CUSTOM_MIE_TABLES_DIR="$MIE_LOCAL/"
        echo "Mie tables: staged $n_copied files ($(du -sh "$MIE_LOCAL" | cut -f1)) from $MIE_SRC to $MIE_LOCAL in $(( $(date +%s) - t_copy )) s"
    else
        echo "WARNING: only $n_copied of $MIE_N_EXPECTED Mie tables copied; libRadtran will read them from /projects" >&2
        rm -rf "$MIE_LOCAL"
    fi
else
    echo "WARNING: no complete Mie-table source found; libRadtran will read the default /projects tables" >&2
fi
# -------------------------------------------------------------------------

# Stagger MATLAB starts across concurrent array tasks so that up to 40
# 40-worker parpools do not all hit MATLAB start-up, preferences and process
# launch at the same instant (task 1 waits 0 s, task 40 waits 780 s).
STARTUP_STAGGER_MOD=${STARTUP_STAGGER_MOD_OVERRIDE:-40}
STARTUP_STAGGER_STEP_SECONDS=${STARTUP_STAGGER_STEP_SECONDS_OVERRIDE:-20}
STARTUP_STAGGER_SECONDS=$(( ((SLURM_ARRAY_TASK_ID - 1) % STARTUP_STAGGER_MOD) * STARTUP_STAGGER_STEP_SECONDS ))
echo "Startup stagger: sleeping ${STARTUP_STAGGER_SECONDS}s before MATLAB launch"
sleep "$STARTUP_STAGGER_SECONDS"

echo "Starting MATLAB job for synthetic clouds ${start_id}..${end_id} at $(date)"

# MATLAB preferences in this task's scratch, to avoid collisions with /home
# and with other concurrent array tasks.
export MATLAB_PREFDIR="$TMPDIR/matlab_prefs/R2024b"
mkdir -p "$MATLAB_PREFDIR"
echo "MATLAB_PREFDIR: $MATLAB_PREFDIR"

# Same MATLAB block as the 20-July-2026 production run: two parpool attempts
# (the second with 8 fewer workers after a 300-420 s pause), then one cloud at
# a time with per-cloud try/catch so a single bad cloud does not end the task.
# Per-cloud start times are logged for wall-time tuning.
# Start-up steps are timestamped ([time] elapsed-s step): the 13-Nov test spent
# ~53 min between MATLAB launch and the first cloud, so each step is timed.
matlab_status=0
time matlab -nodesktop -nodisplay -r "try, t_su = tic; fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'MATLAB up'); p = strsplit(genpath('/projects/anbu8374/Matlab-Research'), pathsep); p = p(~contains(p, '/.') & ~cellfun(@isempty, p)); addpath(strjoin(p, pathsep)); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'path: Matlab-Research (no hidden folders)'); fprintf('MATLAB path entries: %d\n', numel(strsplit(path, pathsep))); addLibRadTran_paths; fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'addLibRadTran_paths'); folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'define_folderPaths_for_HySICS'); for pool_attempt = 1:2, try, fprintf('\n=== PARPOOL START ATTEMPT %d of 2 ===\n', pool_attempt); start_parallel_pool(folder_paths.which_computer); break; catch ME_pool, fprintf(2, '\n[PARPOOL START FAILED attempt %d] %s: %s\n', pool_attempt, ME_pool.identifier, ME_pool.message); old_pool = gcp('nocreate'); if ~isempty(old_pool), delete(old_pool); end; if pool_attempt < 2, slurm_cpus = str2double(getenv('SLURM_CPUS_PER_TASK')); if ~isnan(slurm_cpus) && slurm_cpus > 8, retry_workers = slurm_cpus - 8; setenv('SLURM_CPUS_PER_TASK', num2str(retry_workers)); fprintf(2, 'Retrying parpool with %d workers after a short pause.\n', retry_workers); end; pause(300 + randi(120)); else, rethrow(ME_pool); end; end; end; fprintf('\nWorkers in pool: %d\n', gcp('nocreate').NumWorkers); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'parpool ready'); chunk_settings = struct('n_segments', ${N_SEGMENTS}, 'file_tag', '${FILE_TAG}'); summ = run_synthetic_spectra_chunk('${input_file}', ${start_id}:${end_id}, folder_paths, '${output_dir}', chunk_settings); if summ.n_failed > 0, exit(2); else, exit(0); end; catch ME, fprintf(2, '\n[FATAL MATLAB ERROR] %s: %s\n', ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; exit(1); end"
matlab_status=$?

echo "Finished MATLAB job for synthetic clouds ${start_id}..${end_id} at $(date)"

# Cleanup — INP_OUT, wc and atmmod live under $TMPDIR, so one rm -rf covers
# them along with parpool's JobStorageLocation files.
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

# Prune older scratch dirs (>7 days)
echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true

if [ "$matlab_status" -ne 0 ]; then
    echo "MATLAB exited with status ${matlab_status}; marking this Slurm task as failed."
    exit "$matlab_status"
fi
