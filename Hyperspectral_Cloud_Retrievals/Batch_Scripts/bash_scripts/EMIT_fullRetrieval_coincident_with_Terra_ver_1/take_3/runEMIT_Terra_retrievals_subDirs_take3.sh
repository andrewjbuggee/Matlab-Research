#!/bin/bash

# SLURM Job Array Script to run EMIT droplet profile retrievals
# Each array job processes one subdirectory of coincident EMIT-Terra data
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
#SBATCH --cpus-per-task=25
#SBATCH --mem=90G
#SBATCH --time=23:59:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_Terra_dropProf_subdir_%A_%a
#SBATCH --output=EMIT_Terra_dropProf_subdir_%A_%a.out
#SBATCH --error=EMIT_Terra_dropProf_subdir_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=101-118       # UPDATE: set to the number of subdirectories

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
# runUVSPEC_ver2.m and runMIE.m read this to run uvspec/mie with the GSL
# libraries instead of MATLAB's bundled libstdc++ (avoids GLIBCXX errors).
# Without it, uvspec/mie launch with an empty LD_LIBRARY_PATH and fail.
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b


# ----------------------------------------------------------
# *** DEFINE THE PARENT DIRECTORY CONTAINING SUBDIRECTORIES ***
# *** CANNOT HAVE TRAILING SLASH '/' AT THE END             ***
# ----------------------------------------------------------
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Terra_data/southEast_pacific/likely_stratus"
# ----------------------------------------------------------


# ----------------------------------------------------------
# Define the output directory for this job (can be customized as needed), where
# all the outputs from the MATLAB retrievals will be saved.

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


# Per-task scratch for parpool's JobStorageLocation and the libRadtran INP_OUT/
# wc/atmmod dirs (define_EMIT_dataPath_and_saveFolders puts those under $TMPDIR
# when it is set). Prefer node-local /tmp over shared Lustre to avoid the
# multi-minute "Job Queued" parpool waits seen on the neural-network runs;
# falls back to Lustre if /tmp is unwritable. This replaces the old approach of
# rm-ing the shared ~/.matlab/local_cluster_jobs dir, which could clobber other
# concurrent jobs.
TMPDIR_LOCAL="/tmp/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if mkdir -p "${TMPDIR_LOCAL}" 2>/dev/null && [ -w "${TMPDIR_LOCAL}" ]; then
    export TMPDIR="${TMPDIR_LOCAL}"
    echo "TMPDIR (node-local): ${TMPDIR}"
else
    export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "${TMPDIR}"
    echo "TMPDIR (Lustre fallback): ${TMPDIR}"
fi
echo "TMPDIR headroom:"
df -h "${TMPDIR}"

# Per-task MATLAB preferences dir (node-local) so concurrent array tasks don't
# collide in ~/.matlab or fight over ~/.matlab/local_cluster_jobs.
export MATLAB_PREFDIR="${TMPDIR}/matlab_prefs/R2024b"
mkdir -p "${MATLAB_PREFDIR}"
echo "MATLAB_PREFDIR: ${MATLAB_PREFDIR}"


# Stagger MATLAB startups so concurrent array tasks don't all hit MATLAB+parpool
# launch at the same instant (15 s per task within each wave of 40).
STARTUP_STAGGER_SECONDS=$(( ((SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN) % 40) * 15 ))
echo "Startup stagger: sleeping ${STARTUP_STAGGER_SECONDS}s before MATLAB launch"
sleep "${STARTUP_STAGGER_SECONDS}"


# Run MATLAB
# - define_EMIT_dataPath_and_saveFolders uses SLURM_ARRAY_TASK_ID as the
#   folder_extension_number so each job gets unique INP_OUT, wc, atmmod,
#   and Mie directories
# - folder_paths.coincident_dataFolder is set to the subdirectory name
#   (with trailing slash as expected by the MATLAB code)
echo " "
echo "Starting MATLAB at $(date)"

# try/catch so a failed subdirectory prints a stack trace and returns a non-zero
# exit status (instead of MATLAB exiting 0 on an uncaught error), making the
# task visible as FAILED in sacct for a targeted resubmit.
time matlab -nodesktop -nodisplay -r "\
    try; \
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
            run_dropProf_acpw_retrieval_EMIT_overlap_Terra_ver5(folder_paths, print_status_updates, print_libRadtran_err, \
            plot_figures, save_figures); \
        exit(0); \
    catch ME; \
        fprintf(2, '\n[RETRIEVAL FAILED] %s: %s\n', ME.identifier, ME.message); \
        for k = 1:numel(ME.stack); fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; \
        exit(1); \
    end"
MATLAB_STATUS=$?

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date) (MATLAB status ${MATLAB_STATUS})"
echo "Processed subdirectory: ${CURRENT_SUBDIR}"

# ----------------------------------------------------------
# Cleanup. INP_OUT/wc/atmmod live under $TMPDIR (node-local) so removing TMPDIR
# covers them along with parpool's JobStorageLocation; the Mie dir stays on
# /scratch. A 7-day prune mops up anything left by killed jobs.
# ----------------------------------------------------------
echo "Removing TMPDIR: ${TMPDIR}"
rm -rf "${TMPDIR}"

echo "Removing this task's Mie directory on /scratch..."
rm -rf "/scratch/alpine/${USER}/Mie_Calculations/emit_${SLURM_ARRAY_TASK_ID}" 2>/dev/null || true

echo "Pruning scratch/projects directories older than 7 days..."
find /scratch/alpine/${USER}/                  -maxdepth 1 -name "matlab_tmp_*"  -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/EMIT/             -maxdepth 1 -name "INP_OUT_*"      -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/${USER}/Mie_Calculations/ -maxdepth 1 -name "emit_*"         -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /projects/${USER}/software/libRadtran-2.0.5/data/ -maxdepth 1 -name "wc_*"     -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /projects/${USER}/software/libRadtran-2.0.5/data/ -maxdepth 1 -name "atmmod_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true

# Report a non-zero status to SLURM if the retrieval failed.
exit ${MATLAB_STATUS}
