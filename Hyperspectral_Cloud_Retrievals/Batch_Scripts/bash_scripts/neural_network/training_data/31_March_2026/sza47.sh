#!/bin/bash

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --time=23:59:59   # Request 23 hours and 59 minutes for longer computation
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=110G        # Should be closer to 80% efficiency based on previous runs, but giving some buffer for variability
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza47
#SBATCH --output=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza47_%A_%a.out
#SBATCH --error=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza47_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=501-573%10    # 73 measurements from the ensemble_profiles to process

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

# Define the solar zenith angle for the measurements 
# The range should be linear in cosine space
# mu_samples = linspace(cosd(0), cosd(65), 8);
# acosd(mu_samples) = 0   23.4343   33.3806   41.1882   47.9277   54.0142   59.6619   65.0000
SZA=47.9277

# define the offset for the job array so that the measurement index ranges from 1 to 73 instead of 101 to 173
offset=500

# define the output directory for the results
# ** Needs a trailing slash at the end for the MATLAB script to work properly **
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/dataSet_created_on_31_March_2026/"
mkdir -p $output_dir

# Create unique temp directory for this array task to avoid race conditions
export TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Add a small random delay to prevent simultaneous MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 120))

# Run MATLAB script for the specific measurement index
echo "Starting MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"

# 
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); measurement_idx = ${SLURM_ARRAY_TASK_ID} - ${offset}; hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver3(folder_paths, measurement_idx, ${SZA}, '${output_dir}'); exit"

echo "Finished MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"

# -------------------------------------------------------
# Cleanup: immediately remove this job's MATLAB temp dir.
# This is the JobStorageLocation set in start_parallel_pool.m
# and contains small MATLAB parallel pool bookkeeping files.
# -------------------------------------------------------
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

# Cleanup this job's INP_OUT directory (INP/OUT files already deleted in parfor; removes errMsg.txt + dir)
echo "Removing INP_OUT directory for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/HySICS/INP_OUT_${SLURM_ARRAY_TASK_ID}"

# Cleanup this job's atmmod and wc files older than 7 days (if any remain from previous runs)
echo "Removing atmmod and wc directories for task ${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/atmmod_${SLURM_ARRAY_TASK_ID}"
rm -rf "/scratch/alpine/${USER}/software/libRadtran-2.0.5/data/wc_${SLURM_ARRAY_TASK_ID}"


# Prune scratch directories older than 7 days
echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
