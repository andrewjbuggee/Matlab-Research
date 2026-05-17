#!/bin/bash
#
# Per-chunk Slurm job array for the synthetic-cloud RT training set.
# Each task processes CHUNK_SIZE consecutive clouds (rows of the .nc input)
# inside a SINGLE MATLAB session, reusing the parpool across all clouds in
# the chunk. This requires the matching edit to start_parallel_pool.m so it
# does not tear down a healthy pool on the second-and-later calls.
#
# Cloud index range covered by task t (1-indexed):
#   start = CLOUD_OFFSET + (t-1)*CHUNK_SIZE + 1
#   end   = min(CLOUD_OFFSET + t*CHUNK_SIZE, N_TOTAL)
#
# Submit subsequent batches by bumping CLOUD_OFFSET (and shrinking the
# --array range on the final batch). With CHUNK_SIZE=57 and 1000 tasks per
# batch, six batches cover all 300,001 clouds (last batch: --array=1-248%48
# with CLOUD_OFFSET=285900).
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=82G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=create_meas_synthetic_NN_trainingData_BATCH
#SBATCH --chdir=/projects/anbu8374/slurm_logs/29_April_2026_synthetic
#SBATCH --output=create_meas_synthetic_NN_trainingData_BATCH_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_BATCH_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-900%48

# --- Chunk parameters (edit between batches) -----------------------------
CHUNK_SIZE=57       # clouds per array task
CLOUD_OFFSET=57908    # finished through cloud 57908 in array job 26685405
N_TOTAL=300001      # length of the 'cloud' dimension in the input .nc
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

cd /projects/anbu8374/

module load matlab/R2024b

# Input file (the .nc produced by 07_build_training_inputs.py)
input_file="/scratch/alpine/anbu8374/neural_network_training_data/training_inputs_jointMVN_N300001_L7.nc"

# Output dir — same as synthetic.sh, so this run extends the existing
# 29-April-2026 dataset starting at cloud_id 901. (Trailing slash required.)
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/synthetic_dataSet_created_on_29_April_2026/"
mkdir -p "$output_dir"

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

# Per-task scratch on node-local /tmp instead of Lustre. This is what moves
# parpool's JobStorageLocation off the shared filesystem and avoids the
# multi-minute "Job Queued" parpool waits we saw on tasks 481 and 525 of
# the prior run. /tmp is small but parpool job state is KB-sized.
# Falls back to Lustre if /tmp is unwritable for any reason.
TMPDIR_LOCAL="/tmp/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if mkdir -p "$TMPDIR_LOCAL" 2>/dev/null && [ -w "$TMPDIR_LOCAL" ]; then
    export TMPDIR="$TMPDIR_LOCAL"
    echo "TMPDIR (node-local): $TMPDIR"
else
    export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
    echo "TMPDIR (Lustre fallback): $TMPDIR"
fi

# Report headroom on the chosen TMPDIR so we can spot capacity problems
# in the .out log without having to ssh into the node.
echo "TMPDIR headroom:"
df -h "$TMPDIR"

# Pre-create MATLAB-managed parent paths (Lustre mkdir can fail under load).
# INP_OUT/wc/atmmod now live under $TMPDIR (node-local /tmp), so only the
# Mie parent — which stays on /scratch so tables can be reused across runs —
# needs pre-creating here.
mkdir -p "/scratch/alpine/${USER}/Mie_Calculations/"

# Stagger MATLAB starts across concurrent array tasks.
sleep $((SLURM_ARRAY_TASK_ID % 10))

echo "Starting MATLAB job for synthetic clouds ${start_id}..${end_id} at $(date)"

# Set MATLAB preferences dir to a scratch path to avoid collisions with /home and stop ~/.matlab/local_cluster_jobs from refilling. This also ensures the parpool JobStorageLocation is on the local scratch, which eliminates the multi-minute "Job Queued" waits we saw on tasks 481 and 525 of the prior run. The preferences dir must be writable and persist for the duration of the MATLAB session, but can be shared across tasks since we set a unique JobStorageLocation below.
export MATLAB_PREFDIR=/projects/$USER/.matlab_prefs/R2024b
mkdir -p "$MATLAB_PREFDIR"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); start_parallel_pool(folder_paths.which_computer); input_file = '${input_file}'; output_dir = '${output_dir}'; n_ok = 0; n_fail = 0; for cloud_id = ${start_id}:${end_id}, fprintf('\n=== STARTING cloud %d ===\n', cloud_id); t0 = tic; try, hysics_refl_from_synthetic_NN_inputs(input_file, cloud_id, folder_paths, output_dir); fprintf('=== FINISHED cloud %d in %.1f s ===\n', cloud_id, toc(t0)); n_ok = n_ok + 1; catch ME, fprintf('\n[CLOUD %d FAILED] %s: %s\n', cloud_id, ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf('  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; n_fail = n_fail + 1; end; end; fprintf('\n=== TASK SUMMARY: %d ok, %d failed ===\n', n_ok, n_fail); exit"

echo "Finished MATLAB job for synthetic clouds ${start_id}..${end_id} at $(date)"

# Cleanup — INP_OUT, wc, and atmmod all live under $TMPDIR now, so a single
# rm -rf "$TMPDIR" covers them along with parpool's JobStorageLocation files.
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

# Prune older scratch dirs (>7 days)
echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
