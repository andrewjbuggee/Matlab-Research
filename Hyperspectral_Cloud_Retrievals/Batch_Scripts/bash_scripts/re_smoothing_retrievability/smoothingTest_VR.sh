#!/bin/bash

# Retrievability test of vertical droplet-size structure -- VOCALS-REx
# Each array task processes ONE in-situ profile: it computes the 636-band
# HySICS nadir reflectance for the raw r_e(z) profile plus the smoothed
# versions defined in define_re_smoothing_windows.m, holding tau_c fixed.

# %A: Job ID   %a: Array Task ID
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=reSmoothTest_VR
#SBATCH --output=reSmoothTest_VR_%A_%a.out
#SBATCH --error=reSmoothTest_VR_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=2001-2073%40    # 73 VOCALS-REx profiles (task id = 2000 + profile index)

# Load modules
ml purge
ml gcc/11.2.0
ml netcdf/4.8.1
ml perl/5.36.0
ml texlive/2021

# Set up environment paths
export PATH=/projects/$USER/software/libRadtran-2.0.5/:$PATH
export PATH=/projects/$USER/software/libRadtran-2.0.5/data/:$PATH
export PATH=/projects/$USER/software/libRadtran-2.0.5/bin/:$PATH
export GSL_BIN=/projects/$USER/software/gsl-2.6/bin
export GSL_LIB=/projects/$USER/software/gsl-2.6/lib
export GSL_INC=/projects/$USER/software/gsl-2.6/include
export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH
export INSTALL_DIR=/projects/$USER/software/libRadtran-2.0.5
export PATH=$GSL_BIN:$PATH

# *** Capture the correct LD_LIBRARY_PATH before MATLAB contaminates it ***
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# switch to home directory
cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b

# ---- Experiment parameters ----
# Solar zenith angle for the comparison (single fixed geometry; nadir view).
SZA=30

# Campaign string passed to the MATLAB function
CAMPAIGN=vocalsRex

# Offset so that array tasks 2001-2073 map to profile indices 1-73
offset=2000

# Output directory (trailing slash required by the MATLAB script)
output_dir="/scratch/alpine/anbu8374/re_smoothing_retrievability/dataSet_created_on_23_June_2026/"
mkdir -p $output_dir

# Unique temp directory for this array task to avoid race conditions
export TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Stagger MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 120))

echo "Starting MATLAB job for VOCALS profile (task ${SLURM_ARRAY_TASK_ID}) at $(date)"

# The 6th argument [] tells the function to use define_re_smoothing_windows.m
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); measurement_idx = ${SLURM_ARRAY_TASK_ID} - ${offset}; hysics_refl_smoothingTest_from_insitu(folder_paths, measurement_idx, ${SZA}, '${output_dir}', '${CAMPAIGN}', []); exit"

echo "Finished MATLAB job for VOCALS profile (task ${SLURM_ARRAY_TASK_ID}) at $(date)"

# -------------------------------------------------------
# Cleanup
# -------------------------------------------------------
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

echo "Removing INP_OUT directory for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/HySICS/INP_OUT_${SLURM_ARRAY_TASK_ID}"

echo "Removing atmmod and wc directories for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/atmmod_${SLURM_ARRAY_TASK_ID}"
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/wc_${SLURM_ARRAY_TASK_ID}"

# Prune scratch directories older than 7 days
echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
