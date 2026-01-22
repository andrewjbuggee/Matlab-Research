#!/bin/bash

# Hybrid approach: Process multiple files per job to maximize node utilization

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
#SBATCH --time=1:30:00     # Longer time for multiple files
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=coompute_mieTable_for_libRadtran_gammaDist_wavelength_idx_%A_%a
#SBATCH --output=compute_mieTable_for_libRadtran_gammaDist_wavelength_idx_%A_%a.out
#SBATCH --error=compute_mieTable_for_libRadtran_gammaDist_wavelength_idx_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1001-1002       # 2 jobs × 1 files each = 2 files for subset folder

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



# Define the directory that contains the 2 subdirectories
# ----------------------------------------------------------
# *** MODIFY THIS DIRECTORY BASED ON THE LOCATION OF THE MEASUREMENTS ***
# *** CANNOT HAVE TRAILING SLASH '/' AT THE END         ***
INPUT_DIR="/projects/anbu8374/Matlab-Research/Radiative_Transfer_Physics/mieTables_gamma/netCDF_gammaDist_TEST"
# ----------------------------------------------------------



# Get list of all subdirectories within INPUT_DIR
# ----------------------------------------------------------
# *** NO NEED TO MODIFY THIS SECTION UNLESS SEARCHING CRITERIA DIFFERS ***
mapfile -t ALL_SUBDIRS < <(find "${INPUT_DIR}" -maxdepth 1 -type d ! -path "${INPUT_DIR}" -printf "%f\n" | sort)



# Calculate which directories this job should process
# ----------------------------------------------------------
# *** MODIFY THIS VALUE BASED ON NUMBER OF FILES AND JOBS ***
DIRECTORIES_PER_JOB=1
# ----------------------------------------------------------

START_IDX=$(( (SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN) * DIRECTORIES_PER_JOB ))
END_IDX=$(( START_IDX + DIRECTORIES_PER_JOB - 1 ))

# Total number of subdirectories found
TOTAL_SUBDIRS=${#ALL_SUBDIRS[@]}



# Start of the job
echo " "
echo "Starting MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "Processing directories ${START_IDX} to ${END_IDX}"


# Build the directory path array for this job
DIRECTORY_PATHS=""
for (( i=START_IDX; i<=END_IDX; i++ )); do
    if [ $i -lt ${#ALL_SUBDIRS[@]} ]; then
        CURRENT_SUBDIR=${ALL_SUBDIRS[$i]}
        FULL_PATH="${INPUT_DIR}/${CURRENT_SUBDIR}"
        
        if [[ ! -d "${FULL_PATH}" ]]; then
            echo "WARNING: Directory ${FULL_PATH} does not exist!"
            continue
        fi
        
        # Add directory path to array string (MATLAB cell array format)
        if [ -z "$DIRECTORY_PATHS" ]; then
            DIRECTORY_PATHS="'${FULL_PATH}'"
        else
            DIRECTORY_PATHS="${DIRECTORY_PATHS}, '${FULL_PATH}'"
        fi
        echo "Added to processing list: ${FULL_PATH}"
    fi
done

echo " "
echo "Total directories to process in this job: $(echo $DIRECTORY_PATHS | tr ',' '\n' | wc -l)"

# -------------------------------------------------------------
# Debugging section to check if directory array is defined correctly
echo " "
echo "=== DEBUG INFO ==="
echo "SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Total subdirectories found: ${#ALL_SUBDIRS[@]}"
echo "START_IDX: ${START_IDX}"
echo "END_IDX: ${END_IDX}"
echo "DIRECTORY_PATHS contents: [${DIRECTORY_PATHS}]"

# List the first few subdirectories found:
echo "First few subdirectories found:"
for (( j=0; j<5 && j<${#ALL_SUBDIRS[@]}; j++ )); do
    echo "  [$j]: ${ALL_SUBDIRS[$j]}"
done
echo "=================="

# Check if DIRECTORY_PATHS is empty before running MATLAB
if [ -z "$DIRECTORY_PATHS" ]; then
    echo "ERROR: No directories assigned to this job!"
    echo "This usually means the array indexing is wrong."
    exit 1
fi
# ----------------------------------------------

# Run MATLAB once with all directories for this job
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; precomputeMieTables_slurm({${DIRECTORY_PATHS}}); exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"