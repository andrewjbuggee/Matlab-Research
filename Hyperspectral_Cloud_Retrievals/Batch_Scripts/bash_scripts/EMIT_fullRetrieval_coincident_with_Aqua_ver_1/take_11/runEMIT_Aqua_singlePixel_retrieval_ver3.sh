#!/bin/bash

# SLURM Job Array Script to run EMIT droplet profile retrievals
# MULTIPLE PIXELS PER SLURM TASK
#
# This script is designed to work with per-pixel .mat files created by
# save_overlap_data_perPixel_EMIT_Aqua.m. Each .mat file contains all the
# data needed to run the retrieval for a single EMIT pixel.
#
# The INPUT_DIR should contain .mat files like:
#   overlap_EMIT_pixel_001_2023_9_16_T191118_1.mat
#   overlap_EMIT_pixel_002_2023_9_16_T191118_1.mat
#   etc.
#
# The .mat files are distributed evenly across SLURM array tasks.
# Each task processes ceil(TOTAL_FILES / NUM_JOBS) files (the last task
# may process fewer if the division is not even).
#
# Example: 300 files with --array=1-100 → 3 files per task
# Example: 5 files with --array=1-2   → task 1 gets 3, task 2 gets 2

# %A: Job ID
# %a: Array Task ID

# ----------------------------------------------------------
# *** UPDATE JOB ARRAY RANGE ***
# The array range does NOT need to match the number of .mat files.
# Files are distributed evenly across the array tasks.
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=90G
#SBATCH --time=23:59:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_singlePix_%A_%a
#SBATCH --output=EMIT_singlePix_%A_%a.out
#SBATCH --error=EMIT_singlePix_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1001-1672       # There are 672 files spread across 336 jobs (2 files/job, last job has 2 files)

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


# ----------------------------------------------------------
# *** DEFINE THE DIRECTORY CONTAINING PER-PIXEL .mat FILES ***
# *** CANNOT HAVE TRAILING SLASH '/' AT THE END             ***
# ----------------------------------------------------------
INPUT_DIR="/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver2"
# ----------------------------------------------------------

# ---------------------------------------------------------------
# *** DEFINE THE DIRECTORY WHERE ALL .mat FILES WILL BE SAVED ***
# *** MUST HAVE TRAILING SLASH '/' AT THE END             ***
# ---------------------------------------------------------------
OUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/Droplet_profile_retrievals/take_11/"

mkdir -p "${OUT_DIR}"
# ----------------------------------------------------------



# Get sorted list of .mat files (full paths)
mapfile -t ALL_MAT_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f | sort)

TOTAL_FILES=${#ALL_MAT_FILES[@]}

if [ ${TOTAL_FILES} -eq 0 ]; then
    echo "ERROR: No .mat files found in ${INPUT_DIR}"
    exit 1
fi

# Calculate the number of jobs in the array
SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN:-1001}  # Default to 1 if unset
NUM_JOBS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))

# Calculate files per job (ceiling division)
FILES_PER_JOB=$(( (TOTAL_FILES + NUM_JOBS - 1) / NUM_JOBS ))

# Calculate the start and end indices for this task (0-based)
JOB_IDX=$(( SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN ))
START_IDX=$(( JOB_IDX * FILES_PER_JOB ))
END_IDX=$(( START_IDX + FILES_PER_JOB - 1 ))

# Clamp END_IDX to the last file
if [ ${END_IDX} -ge ${TOTAL_FILES} ]; then
    END_IDX=$(( TOTAL_FILES - 1 ))
fi

# If START_IDX is beyond the file list, this task has nothing to do
if [ ${START_IDX} -ge ${TOTAL_FILES} ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} has no files to process (${TOTAL_FILES} files, ${NUM_JOBS} jobs, ${FILES_PER_JOB} files/job)"
    echo "Exiting gracefully."
    exit 0
fi

N_FILES_THIS_JOB=$(( END_IDX - START_IDX + 1 ))


# Start of the job
echo " "
echo "Starting MATLAB multi-pixel retrieval at $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "Total .mat files: ${TOTAL_FILES}, Jobs: ${NUM_JOBS}, Files/job: ${FILES_PER_JOB}"
echo "This task processes files ${START_IDX} to ${END_IDX} (${N_FILES_THIS_JOB} files)"


# -------------------------------------------------------------
# Debugging section
echo " "
echo "=== DEBUG INFO ==="
echo "SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Total .mat files found: ${TOTAL_FILES}"
echo "NUM_JOBS: ${NUM_JOBS}"
echo "FILES_PER_JOB: ${FILES_PER_JOB}"
echo "JOB_IDX: ${JOB_IDX}"
echo "START_IDX: ${START_IDX}"
echo "END_IDX: ${END_IDX}"


# ----------------------------------------------------------
# This code below is repeated in each loop iteration to ensure a clean MATLAB environment for each file
# Give each SLURM task its own isolated MATLAB preferences and job directory
MATLAB_TASK_DIR="/scratch/alpine/${USER}/matlab_jobs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${MATLAB_TASK_DIR}"
export MATLAB_TASK_DIR
echo "MATLAB_TASK_DIR: ${MATLAB_TASK_DIR}"
# ----------------------------------------------------------
# ----------------------------------------------------------


# List first 10 .mat files found:
echo " "
echo "First 10 .mat files found:"
for (( j=0; j<10 && j<${#ALL_MAT_FILES[@]}; j++ )); do
    echo "  [$j]: ${ALL_MAT_FILES[$j]}"
done

# List the files assigned to this task:
echo " "
echo "Files assigned to this task:"
for (( j=START_IDX; j<=END_IDX; j++ )); do
    echo "  [$j]: ${ALL_MAT_FILES[$j]}"
done
echo "=================="



# Add a small random delay to prevent simultaneous MATLAB startups
# This is especially important if many tasks start at the same time, as it can prevent overwhelming the filesystem and reduce contention for resources.
sleep $((SLURM_ARRAY_TASK_ID % 60))


# Loop over all .mat files assigned to this task
for (( FILE_IDX=START_IDX; FILE_IDX<=END_IDX; FILE_IDX++ )); do

    CURRENT_MAT_FILE=${ALL_MAT_FILES[$FILE_IDX]}

    # Verify the file exists
    if [[ ! -f "${CURRENT_MAT_FILE}" ]]; then
        echo "WARNING: .mat file ${CURRENT_MAT_FILE} does not exist! Skipping."
        continue
    fi

    # Compute a unique folder_extension_number for each file within this task
    # This ensures each call gets unique INP_OUT, wc, atmmod, and Mie directories
    # Formula: TASK_ID * 1000 + local file index (0-based within this task)
    LOCAL_FILE_IDX=$(( FILE_IDX - START_IDX ))
    FOLDER_EXT_NUM=$(( SLURM_ARRAY_TASK_ID * 1000 + LOCAL_FILE_IDX ))

    echo " "
    echo "=============================================="
    echo "Processing file $(( LOCAL_FILE_IDX + 1 )) of ${N_FILES_THIS_JOB}"
    echo "File index: ${FILE_IDX}"
    echo "Mat file: ${CURRENT_MAT_FILE}"
    echo "folder_extension_number: ${FOLDER_EXT_NUM}"
    echo "Starting MATLAB at $(date)"
    echo "=============================================="

    time matlab -nodesktop -nodisplay -r "\
        addpath(genpath('/projects/anbu8374/Matlab-Research')); \
        addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); \
        addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); \
        clear variables; \
        addLibRadTran_paths; \
        print_status_updates = true; \
        print_libRadtran_err = false; \
        mat_file_path = '${CURRENT_MAT_FILE}'; \
        folder_extension_number = ${FOLDER_EXT_NUM}; \
        output_dir = '${OUT_DIR}'; \
        [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = \
            run_retrieval_singlePixel_EMIT_Aqua(mat_file_path, folder_extension_number, \
            print_status_updates, print_libRadtran_err, output_dir); \
        exit"

    echo " "
    echo "Finished processing: ${CURRENT_MAT_FILE} at $(date)"

done

# Clean MATLAB temp directories
echo "Cleaning MATLAB temp directories for task ${SLURM_ARRAY_TASK_ID}"
# Clean up this task's unique MATLAB job directory after completion (at end of script)
rm -rf "${MATLAB_TASK_DIR}"

echo " "
echo "=============================================="
echo "Task ${SLURM_ARRAY_TASK_ID} finished all ${N_FILES_THIS_JOB} files at $(date)"
echo "=============================================="


