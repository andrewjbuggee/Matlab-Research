#!/bin/bash
#
# TEST run of the synthetic-cloud RT pipeline on an Alpine testing node,
# for the 12-Nov-2026 training inputs (K_re = 3, K_lwc = 5; libRadtran phi0
# azimuth convention; mixed-layer water-cloud files).
#
# One array task runs a short, explicit list of cloud indices in a single
# MATLAB session (same driver and parpool logic as the production array).
# The list is in priority order, so whatever finishes inside the 1-hour
# testing limit is the most informative:
#
#     553    typical cloud      (tau_c = 7.8, the median)
#   15022    edge case          LWC = 0 at the cloud-top level
#   27045    edge case          cloud-top LWC = 8.97 g/m^3 (largest in the set)
#    1105    thinnest cloud     tau_c = 3.00
#   21014    thickest cloud     tau_c = 39.97
#
# Indices are 1-based rows of the 'cloud' dimension (the idx the MATLAB driver
# takes). Override the list without editing:  CLOUD_IDS_OVERRIDE="553 15022" sbatch ...
#
# BEFORE SUBMITTING (on Alpine):
#   1. Code: git pull in /projects/anbu8374/Matlab-Research so it has commit
#      6a30579 (write_wc_file_mixed_layers.m). The pre-flight block below
#      prints the commit and stops if that file is missing.
#   2. Input: copy the new .nc from the laptop, e.g.
#        scp lasp-CU-paper-3/training_inputs/training_inputs_jointMVN_N100000_L7_2026-11-12.nc \
#            anbu8374@login.rc.colorado.edu:/scratch/alpine/anbu8374/neural_network_training_data/
#   3. Log directory: Slurm cannot create the --chdir directory, so
#        mkdir -p /projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch_test
#   4. Submit from a login node:  sbatch test_synthetic_updates.sh
#
# Resources (1 h, 10 CPUs, 50 GB) match the earlier test scripts in
# 29_April_2026_synthetic/. Those used --partition=atesting --qos=testing, but
# by Nov 2026 Alpine's `sinfo -s` lists acpu (default CPU partition, 1-day
# limit) and no atesting or amilan, so this test runs on acpu. acpu only
# accepts the QoS values cpu-normal, cpu-long and cpu-testing (sbatch error
# message); testing jobs use cpu-testing. If sbatch rejects the 10-CPU / 50 GB /
# 1 h request under cpu-testing, check its limits with
#   sacctmgr show qos cpu-testing format=name,maxwall,maxtrespu%40
# and shrink to fit, or switch to cpu-normal.
# Submit from a login node after `module load slurm/alpine`.
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=acpu        # Alpine CPU partition (Nov 2026 sinfo lists no atesting or amilan)
#SBATCH --qos=cpu-testing      # acpu accepts cpu-normal / cpu-long / cpu-testing (sbatch error, 12 Nov 2026)
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=create_meas_synthetic_NN_trainingData_TEST_12Nov
#SBATCH --chdir=/projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch_test
#SBATCH --output=create_meas_synthetic_NN_trainingData_TEST_12Nov_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_TEST_12Nov_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1

# --- Test parameters -----------------------------------------------------
CLOUD_IDS=${CLOUD_IDS_OVERRIDE:-"553 15022 27045 1105 21014"}
N_TOTAL=${N_TOTAL_OVERRIDE:-100000}      # length of the 'cloud' dimension in the input .nc
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

# Test outputs go to their own folder so they never mix with a production set.
# (Trailing slash required by the MATLAB driver.)
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/synthetic_dataSet_TEST_12_Nov_2026/"
mkdir -p "$output_dir"

# --- Pre-flight checks: fail in seconds, not after MATLAB starts -----------
code_dir="/projects/anbu8374/Matlab-Research"
echo "Matlab-Research commit: $(git -C "$code_dir" log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo 'unknown (not a git checkout?)')"
if [ -z "$(find "$code_dir" -name write_wc_file_mixed_layers.m -print -quit)" ]; then
    echo "ERROR: write_wc_file_mixed_layers.m not found under $code_dir. git pull first." >&2
    exit 1
fi
if [ ! -f "$input_file" ]; then
    echo "ERROR: input file not found: $input_file (scp it from the laptop first)." >&2
    exit 1
fi
for id in $CLOUD_IDS; do
    if [ "$id" -lt 1 ] || [ "$id" -gt "$N_TOTAL" ]; then
        echo "ERROR: cloud id $id is outside 1..$N_TOTAL" >&2
        exit 1
    fi
done
# MATLAB vector literal, e.g. "553 15022 27045" -> "[553 15022 27045]"
cloud_id_vec="[${CLOUD_IDS}]"
echo "Task ${SLURM_ARRAY_TASK_ID}: test clouds ${cloud_id_vec}"
echo "Input : $input_file"
echo "Output: $output_dir"
# -------------------------------------------------------------------------

# Per-task scratch on node-local /tmp instead of Lustre (keeps parpool's
# JobStorageLocation off the shared filesystem). Falls back to Lustre if /tmp
# is unwritable.
TMPDIR_LOCAL="/tmp/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if mkdir -p "$TMPDIR_LOCAL" 2>/dev/null && [ -w "$TMPDIR_LOCAL" ]; then
    export TMPDIR="$TMPDIR_LOCAL"
    echo "TMPDIR (node-local): $TMPDIR"
else
    export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
    echo "TMPDIR (Lustre fallback): $TMPDIR"
fi
echo "TMPDIR headroom:"
df -h "$TMPDIR"

# Mie parent stays on /scratch so tables can be reused across runs.
mkdir -p "/scratch/alpine/${USER}/Mie_Calculations/"

# (No start-up stagger: a single test task has nothing to collide with.)

echo "Starting MATLAB test job at $(date)"

# MATLAB preferences in this task's scratch to avoid collisions with /home.
export MATLAB_PREFDIR="$TMPDIR/matlab_prefs/R2024b"
mkdir -p "$MATLAB_PREFDIR"
echo "MATLAB_PREFDIR: $MATLAB_PREFDIR"

# Same MATLAB block as the production array, with two changes for the test:
#   * loops over the explicit cloud_id_vec instead of start_id:end_id;
#   * the parpool retry pauses 60 s (production waits 300-420 s, which would
#     eat most of the 1-hour testing limit).
matlab_status=0
time matlab -nodesktop -nodisplay -r "try, addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); for pool_attempt = 1:2, try, fprintf('\n=== PARPOOL START ATTEMPT %d of 2 ===\n', pool_attempt); start_parallel_pool(folder_paths.which_computer); break; catch ME_pool, fprintf(2, '\n[PARPOOL START FAILED attempt %d] %s: %s\n', pool_attempt, ME_pool.identifier, ME_pool.message); old_pool = gcp('nocreate'); if ~isempty(old_pool), delete(old_pool); end; if pool_attempt < 2, pause(60); else, rethrow(ME_pool); end; end; end; fprintf('\nWorkers in pool: %d\n', gcp('nocreate').NumWorkers); input_file = '${input_file}'; output_dir = '${output_dir}'; n_ok = 0; n_fail = 0; for cloud_id = ${cloud_id_vec}, fprintf('\n=== STARTING cloud %d at %s ===\n', cloud_id, char(datetime('now'))); t0 = tic; try, hysics_refl_from_synthetic_NN_inputs(input_file, cloud_id, folder_paths, output_dir); fprintf('=== FINISHED cloud %d in %.1f s ===\n', cloud_id, toc(t0)); n_ok = n_ok + 1; catch ME, fprintf('\n[CLOUD %d FAILED after %.1f s] %s: %s\n', cloud_id, toc(t0), ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf('  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; n_fail = n_fail + 1; end; end; fprintf('\n=== TASK SUMMARY: %d ok, %d failed ===\n', n_ok, n_fail); if n_fail > 0, exit(2); else, exit(0); end; catch ME, fprintf(2, '\n[FATAL MATLAB ERROR] %s: %s\n', ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; exit(1); end"
matlab_status=$?

echo "Finished MATLAB test job at $(date)"

# What the test produced (newest first), so the .out log shows it directly.
echo "Files written to $output_dir:"
ls -lt "$output_dir" | head -20

# Cleanup — INP_OUT, wc and atmmod live under $TMPDIR.
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

# (The production script's 7-day scratch pruning is left out of the test.)

if [ "$matlab_status" -ne 0 ]; then
    echo "MATLAB exited with status ${matlab_status}; marking this Slurm task as failed."
    exit "$matlab_status"
fi
