#!/bin/bash

# Hybrid approach: Process multiple files per job to maximize node utilization
# In the directory below, there are 20 files
# With 10 jobs and 20 files, each job will process 2 files

# SLURM Job Array Script to run MATLAB retrievals on multiple files in parallel
# This will run the same analysis on multiple files within a specified directory

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID

# ----------------------------------------------------------
# *** UPDATE JOB NAME, OUTPUT, AND ERROR FILE BASED ON SIM ***
# *** UPDATE JOB ARRAY RANGE BASED ON NUMBER OF FILES  ***
# ----------------------------------------------------------
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --time=23:59:00     # Longer time for multiple files
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=full_retrieval_EMIT_log_newCov_allBands_pix1_rev_2_%A_%a
#SBATCH --output=full_retrieval_EMIT_log_newCov_allBands_pix1_rev_2_%A_%a.out
#SBATCH --error=full_retrieval_EMIT_log_newCov_allBands_pix1_rev_2_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=201       # 1 jobs × 1 files each = 1 files for subset folder

# ** Is there a way for the --array to read the number of files in a directory automatically? **

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

cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b

# Add a small random delay to prevent simultaneous MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 10))

# Run MATLAB script for the specific measurement index
echo "Starting MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"

# 
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); retrieve_drop_prof_H2O_EMIT_overlap_Aqua_ver5; exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"