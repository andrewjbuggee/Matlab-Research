#!/bin/bash

# SLURM Job Array for TBLUT retrieval analysis
# This will run the same analysis on multiple files within a specified directory

# Becayse this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=28G
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test_retrieval_array_%A_%a
#SBATCH --output=test_retrieval_array_%A_%a.out
#SBATCH --error=test_retrieval_array_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-4  # Adjust this range based on number of files/folders

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

# Change to your working directory
cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b

# Start of the job
echo "Starting MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"

# Define the directory containing your input files
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/test_directory"


# Method 2: Create array from actual files in your directory
# Get list of files (adjust the pattern as needed - e.g., *.mat, *.txt, etc.)
mapfile -t FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f -printf "%f\n" | sort)
# Alternative if you want all files: mapfile -t FILES < <(ls "${INPUT_DIR}")

# Get the current file for this array task
CURRENT_FILE=${FILES[$((SLURM_ARRAY_TASK_ID-1))]}
CURRENT_FILE_PATH="${INPUT_DIR}/${CURRENT_FILE}"

echo "Processing file: ${CURRENT_FILE}"
echo "Full path: ${CURRENT_FILE_PATH}"

# Check if file exists
if [[ ! -f "${CURRENT_FILE_PATH}" ]]; then
    echo "ERROR: File ${CURRENT_FILE_PATH} does not exist!"
    exit 1
fi

# Run your MATLAB script with the current file as an argument
# Each task processes a different file based on the array index
time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); 
clear variables; addLibRadTran_paths; 
current_file='${CURRENT_FILE}';
folder_paths = define_folderPaths_for_HySICS('${SLURM_ARRAY_TASK_ID}'); which_computer = folder_paths.which_computer; 
print_status_updates = true; print_libRadtran_err = true; test_retrieval_HySICS_jobArray('${CURRENT_FILE}'); exit"

echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"