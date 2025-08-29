#!/bin/bash

# Hybrid approach: Process multiple files per job to maximize node utilization
# In the directory below, there are 792 files
# With 36 jobs and 792 files, each job will process 22 files

# SLURM Job Array Script to run MATLAB retrievals on multiple files in parallel
# This will run the same analysis on multiple files within a specified directory

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28G
#SBATCH --time=23:59:59     # Longer time for multiple files
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=full_retrieval_hysics_rTop_10_%A_%a
#SBATCH --output=full_retrieval_hysics_rTop_10_%A_%a.out
#SBATCH --error=full_retrieval_hysics_rTop_10_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-36        # 36 jobs × 22 files each = 792 files

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

# Define the directory containing your input files
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91"

# Get list of all files
mapfile -t ALL_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f -printf "%f\n" | sort)

# Calculate which files this job should process
FILES_PER_JOB=22
START_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) * FILES_PER_JOB ))
END_IDX=$(( START_IDX + FILES_PER_JOB - 1 ))

# Start of the job
echo " "
echo "Starting MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "Processing files ${START_IDX} to ${END_IDX}"

# Build the filename array for this job
FILE_ARRAY=""
for (( i=START_IDX; i<=END_IDX; i++ )); do
    if [ $i -lt ${#ALL_FILES[@]} ]; then
        CURRENT_FILE=${ALL_FILES[$i]}
        if [[ ! -f "${INPUT_DIR}/${CURRENT_FILE}" ]]; then
            echo "WARNING: File ${INPUT_DIR}/${CURRENT_FILE} does not exist!"
            continue
        fi
        
        # Add filename to array string (MATLAB cell array format)
        if [ -z "$FILE_ARRAY" ]; then
            FILE_ARRAY="'${CURRENT_FILE}'"
        else
            FILE_ARRAY="${FILE_ARRAY}, '${CURRENT_FILE}'"
        fi
        echo "Added to processing list: ${CURRENT_FILE}"
    fi
done

echo " "
echo "Total files to process in this job: $(echo $FILE_ARRAY | tr ',' '\n' | wc -l)"

# Run MATLAB once with all files for this job
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS('${SLURM_ARRAY_TASK_ID}'); folder_paths.HySICS_simulated_spectra = '${INPUT_DIR}/'; print_status_updates = true; print_libRadtran_err = true; file_list = {${FILE_ARRAY}}; [tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver3(file_list, folder_paths, print_status_updates, print_libRadtran_err); exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"