#!/bin/bash

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID
# ----------------------------------------------------------

#SBATCH --nodes=1
#SBATCH --time=00:30:00   # Adjusted time for single measurement
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=create_meas_pt3_percent
#SBATCH --output=measurements_pt3_percent_%A_%a.out
#SBATCH --error=create_meas_pt3_percent_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-73    # 73 measurements from the ensemble_profiles to process

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

# switch to home directory
cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b

# Create unique temp directory for this array task to avoid race conditions
export TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Add a small random delay to prevent simultaneous MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 10))

# Run MATLAB script for the specific measurement index
echo "Starting MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"

# 
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); measurement_idx = ${SLURM_ARRAY_TASK_ID}; hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx); exit"

echo "Finished MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"
