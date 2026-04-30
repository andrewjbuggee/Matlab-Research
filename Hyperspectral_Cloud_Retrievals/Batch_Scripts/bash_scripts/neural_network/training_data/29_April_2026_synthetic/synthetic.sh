#!/bin/bash
#
# Per-cloud Slurm job array for the synthetic-cloud RT training set.
# Each task processes one row (cloud) of the .nc input file.
#
# Update --array=1-N to match the number of clouds in input_file.
# (cloud index = SLURM_ARRAY_TASK_ID — no offset needed since the .nc
# is 1-indexed in MATLAB.)
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=31G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=create_meas_synthetic_NN_trainingData
#SBATCH --output=create_meas_synthetic_NN_trainingData_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-900%50     # set to N_clouds in the .nc

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

# Output directory (trailing slash required)
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/synthetic_dataSet_created_on_29_April_2026/"
mkdir -p "$output_dir"

# Per-task scratch dirs (mirror the VOCALS bash patterns)
export TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

sleep $((SLURM_ARRAY_TASK_ID % 10))

echo "Starting MATLAB job for synthetic cloud ${SLURM_ARRAY_TASK_ID} at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); hysics_refl_from_synthetic_NN_inputs('${input_file}', ${SLURM_ARRAY_TASK_ID}, folder_paths, '${output_dir}'); exit"

echo "Finished MATLAB job for synthetic cloud ${SLURM_ARRAY_TASK_ID} at $(date)"

# Cleanup (same pattern as sza0.sh)
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"
echo "Removing INP_OUT directory for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/HySICS/INP_OUT_${SLURM_ARRAY_TASK_ID}"
echo "Removing atmmod and wc directories for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/atmmod_${SLURM_ARRAY_TASK_ID}"
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/wc_${SLURM_ARRAY_TASK_ID}"

# Prune older scratch dirs (>7 days)
echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
