#!/bin/bash

# SLURM Job Array Script to run EMIT droplet profile retrievals
# Each array job processes one subdirectory of coincident EMIT-Aqua data
#
# The input directory contains subdirectories like:
#   2023_9_16_T191118_1/
#   2024_5_17_T183930_1/
#   etc.
#
# Each subdirectory is assigned to a separate job via SLURM job arrays.
# The subdirectory name is passed to MATLAB as folder_paths.coincident_dataFolder

# %A: Job ID
# %a: Array Task ID

# ----------------------------------------------------------
# *** UPDATE JOB ARRAY RANGE TO MATCH NUMBER OF SUBDIRS ***
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G
#SBATCH --time=23:59:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_dropProf_subdir_%A_%a
#SBATCH --output=EMIT_dropProf_subdir_%A_%a.out
#SBATCH --error=EMIT_dropProf_subdir_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-10       # UPDATE: set to the number of subdirectories

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
# *** DEFINE THE PARENT DIRECTORY CONTAINING SUBDIRECTORIES ***
# *** CANNOT HAVE TRAILING SLASH '/' AT THE END             ***
# ----------------------------------------------------------
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific"
# ----------------------------------------------------------


# Get sorted list of subdirectories (names only, no trailing slash)
mapfile -t ALL_SUBDIRS < <(find "${INPUT_DIR}" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | sort)

TOTAL_SUBDIRS=${#ALL_SUBDIRS[@]}

# Calculate the index for this job (0-based)
SUBDIR_IDX=$(( SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN ))

# Check bounds
if [ ${SUBDIR_IDX} -ge ${TOTAL_SUBDIRS} ]; then
    echo "ERROR: Array task ${SLURM_ARRAY_TASK_ID} exceeds number of subdirectories (${TOTAL_SUBDIRS})"
    exit 1
fi

# Get the subdirectory for this job
CURRENT_SUBDIR=${ALL_SUBDIRS[$SUBDIR_IDX]}

# Verify the subdirectory exists
if [[ ! -d "${INPUT_DIR}/${CURRENT_SUBDIR}" ]]; then
    echo "ERROR: Subdirectory ${INPUT_DIR}/${CURRENT_SUBDIR} does not exist!"
    exit 1
fi


# Start of the job
echo " "
echo "Starting MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "Processing subdirectory: ${CURRENT_SUBDIR}"


# -------------------------------------------------------------
# Debugging section
echo " "
echo "=== DEBUG INFO ==="
echo "SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Total subdirectories found: ${TOTAL_SUBDIRS}"
echo "SUBDIR_IDX: ${SUBDIR_IDX}"
echo "CURRENT_SUBDIR: ${CURRENT_SUBDIR}"

# List all subdirectories found:
echo " "
echo "All subdirectories found:"
for (( j=0; j<${#ALL_SUBDIRS[@]}; j++ )); do
    echo "  [$j]: ${ALL_SUBDIRS[$j]}"
done
echo "=================="


# Clean MATLAB temp directories
echo "Cleaning MATLAB temp directories for task ${SLURM_ARRAY_TASK_ID}"
rm -rf ~/.matlab/local_cluster_jobs/R2024b/Job*
rm -rf /tmp/mathworks_${USER}_*


# Add a small random delay to prevent simultaneous MATLAB startups
sleep $((SLURM_ARRAY_TASK_ID % 10))


# Run MATLAB
# - define_EMIT_dataPath_and_saveFolders uses SLURM_ARRAY_TASK_ID as the
#   folder_extension_number so each job gets unique INP_OUT, wc, atmmod,
#   and Mie directories
# - folder_paths.coincident_dataFolder is set to the subdirectory name
#   (with trailing slash as expected by the MATLAB code)
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "\
    addpath(genpath('/projects/anbu8374/Matlab-Research')); \
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); \
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); \
    clear variables; \
    addLibRadTran_paths; \
    print_status_updates = true; \
    print_libRadtran_err = false; \
    plot_figures = false; \
    save_figures = false; \
    folder_paths = define_EMIT_dataPath_and_saveFolders(${SLURM_ARRAY_TASK_ID}); \
    folder_paths.coincident_dataPath = '${INPUT_DIR}/'; \
    folder_paths.coincident_dataFolder = '${CURRENT_SUBDIR}/'; \
    [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = \
        run_dropProf_acpw_retrieval_EMIT_overlap_Aqua_ver5(folder_paths, print_status_updates, print_libRadtran_err, \
        plot_figures, save_figures); \
    exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Processed subdirectory: ${CURRENT_SUBDIR}"
