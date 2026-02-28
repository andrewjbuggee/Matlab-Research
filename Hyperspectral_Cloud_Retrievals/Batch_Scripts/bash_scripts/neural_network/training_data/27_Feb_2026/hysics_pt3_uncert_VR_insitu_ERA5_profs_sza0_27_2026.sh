#!/bin/bash

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --time=23:59:59   # Request 80 hours for longer computation
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza0
#SBATCH --output=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza0_%A_%a.out
#SBATCH --error=create_meas_pt3_percent_VR_insitu_ERA5_trainingData_sza0_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=101-173    # 73 measurements from the ensemble_profiles to process

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

# Define the solar zenith angle for the measurements (0 degrees for this case)
SZA=0

# define the offset for the job array so that the measurement index ranges from 1 to 73 instead of 101 to 173
offset=100

# define the output directory for the results
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/dataSet_created_on_27_Feb_2026"

# Create unique temp directory for this array task to avoid race conditions
export TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Add a small random delay to prevent simultaneous MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 10))

# Run MATLAB script for the specific measurement index
echo "Starting MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"

# 
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); measurement_idx = ${SLURM_ARRAY_TASK_ID} - ${offset}; hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver2(folder_paths, measurement_idx, ${SZA}, '${output_dir}'); exit"

echo "Finished MATLAB job for measurement ${SLURM_ARRAY_TASK_ID} at $(date)"
