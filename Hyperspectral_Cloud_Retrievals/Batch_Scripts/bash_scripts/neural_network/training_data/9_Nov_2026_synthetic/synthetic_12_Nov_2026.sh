#!/bin/bash
#
# PRODUCTION per-chunk Slurm job array for the 12-Nov-2026 synthetic-cloud RT
# training set (100,000 clouds; K_re = 3, K_lwc = 5; libRadtran phi0 azimuth
# convention; mixed-layer water-cloud files). Tested with
# test_synthetic_updates.sh on acpu / cpu-testing.
#
# Each task processes CHUNK_SIZE consecutive clouds (rows of the .nc input)
# inside a SINGLE MATLAB session, reusing the parpool across all clouds in
# the chunk (requires the start_parallel_pool.m that keeps a healthy pool).
#
# Cloud index range covered by task t (1-indexed):
#   start = CLOUD_OFFSET + (t-1)*CHUNK_SIZE + 1
#   end   = min(CLOUD_OFFSET + t*CHUNK_SIZE, N_TOTAL)
#
# Limits and throttling are those of the 20-July-2026 run
# (29_April_2026_synthetic/synthetic_20_July_2026.sh): 40 CPUs, 82 GB, 24 h,
# 57 clouds per task, 900 tasks per batch, at most 40 tasks running at once.
# With CHUNK_SIZE = 57, 100,000 clouds need ceil(100000/57) = 1755 tasks,
# i.e. two batches:
#
#   batch 1:  CLOUD_OFFSET=0       --array=1-900%40   clouds      1 .. 51300
#   batch 2:  CLOUD_OFFSET=51300   --array=1-855%40   clouds  51301 .. 100000
#             (task 855 is clamped at N_TOTAL and runs 22 clouds)
#
# Submit batch 2 without editing the file:
#   CLOUD_OFFSET_OVERRIDE=51300 sbatch --array=1-855%40 synthetic_12_Nov_2026.sh
# (command-line --array overrides the #SBATCH line below).
#
# CHUNK SIZE vs. WALL TIME: 57 clouds per 24 h was tuned on amilan nodes.
# Alpine's CPU partition is now acpu (Nov 2026), on different hardware. Before
# batch 1, take the per-cloud time from the test log ("=== FINISHED cloud N in
# X s ===", on 10 CPUs) and check that 57 clouds fit in ~20 h on 40 CPUs; if
# not, lower CHUNK_SIZE and redo the batch arithmetic above. A task that hits
# the 24-h limit loses its unfinished clouds (each finished cloud is saved).
#
# BEFORE SUBMITTING (on Alpine):
#   1. module load slurm/alpine      (on login-ci4 this was needed for sinfo to
#                                     list acpu; the CURC batch-jobs page does not
#                                     mention it, so it may be harmless if absent)
#   2. git pull in /projects/anbu8374/Matlab-Research   (needs commit 6a30579)
#   3. The input .nc is on /scratch (the test job already used it)
#   4. mkdir -p /projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch1
#      (Slurm cannot create the --chdir directory; use ..._batch2 for batch 2
#       and pass --chdir on the sbatch command line)
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
#SBATCH --time=24:00:00
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --mem=82G               # below 40 cores x 3.8 GB/core = 152 GB, so no extra cores are charged
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=create_meas_synthetic_NN_trainingData_12Nov
#SBATCH --chdir=/projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch1
#SBATCH --output=create_meas_synthetic_NN_trainingData_12Nov_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_12Nov_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-900%40

# --- Chunk parameters (edit or override between batches) -----------------
CHUNK_SIZE=${CHUNK_SIZE_OVERRIDE:-57}         # clouds per array task
CLOUD_OFFSET=${CLOUD_OFFSET_OVERRIDE:-0}      # batch 1: 0; batch 2: 51300
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
input_file="/scratch/alpine/anbu8374/neural_network_training_data/training_inputs_jointMVN_N100000_L7_2026-11-12.nc"

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

# Mie parent stays on /scratch so tables can be reused across runs.
mkdir -p "/scratch/alpine/${USER}/Mie_Calculations/"

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
matlab_status=0
time matlab -nodesktop -nodisplay -r "try, addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); for pool_attempt = 1:2, try, fprintf('\n=== PARPOOL START ATTEMPT %d of 2 ===\n', pool_attempt); start_parallel_pool(folder_paths.which_computer); break; catch ME_pool, fprintf(2, '\n[PARPOOL START FAILED attempt %d] %s: %s\n', pool_attempt, ME_pool.identifier, ME_pool.message); old_pool = gcp('nocreate'); if ~isempty(old_pool), delete(old_pool); end; if pool_attempt < 2, slurm_cpus = str2double(getenv('SLURM_CPUS_PER_TASK')); if ~isnan(slurm_cpus) && slurm_cpus > 8, retry_workers = slurm_cpus - 8; setenv('SLURM_CPUS_PER_TASK', num2str(retry_workers)); fprintf(2, 'Retrying parpool with %d workers after a short pause.\n', retry_workers); end; pause(300 + randi(120)); else, rethrow(ME_pool); end; end; end; fprintf('\nWorkers in pool: %d\n', gcp('nocreate').NumWorkers); input_file = '${input_file}'; output_dir = '${output_dir}'; n_ok = 0; n_fail = 0; for cloud_id = ${start_id}:${end_id}, fprintf('\n=== STARTING cloud %d at %s ===\n', cloud_id, char(datetime('now'))); t0 = tic; try, hysics_refl_from_synthetic_NN_inputs(input_file, cloud_id, folder_paths, output_dir); fprintf('=== FINISHED cloud %d in %.1f s ===\n', cloud_id, toc(t0)); n_ok = n_ok + 1; catch ME, fprintf('\n[CLOUD %d FAILED after %.1f s] %s: %s\n', cloud_id, toc(t0), ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf('  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; n_fail = n_fail + 1; end; end; fprintf('\n=== TASK SUMMARY: %d ok, %d failed ===\n', n_ok, n_fail); if n_fail > 0, exit(2); else, exit(0); end; catch ME, fprintf(2, '\n[FATAL MATLAB ERROR] %s: %s\n', ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; exit(1); end"
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
