#!/bin/bash

# ======================================================================
# SINGLE-PIXEL SMOKE TEST for the take_13 / ver5 EMIT-Aqua retrieval.
#
# Processes EXACTLY ONE per-pixel .mat file, end-to-end, so you can
# confirm the updated pipeline works before launching all 672 tasks.
# It validates the things that changed since take_12:
#   - PRE_MATLAB_LD_LIBRARY_PATH is exported, so uvspec AND mie find libgsl
#     (runUVSPEC_ver2.m / runMIE.m run them via `env LD_LIBRARY_PATH=...`)
#   - the parallel pool starts (SLURM_CPUS_PER_TASK workers, TMPDIR storage)
#   - the full Gauss-Newton retrieval runs and writes its output .mat
#   - end-of-job cleanup (scratch + /projects wc_*/atmmod_*) runs
#
# Submit with:  sbatch TEST_runEMIT_Aqua_singlePixel_retrieval_ver5.sh
#
# Resources are identical to the production script so the test exercises
# the same configuration. The only differences are: a single-element
# array, a distinct job-name/output, ONE file processed, and a separate
# test output directory.
# ======================================================================

# %A: Job ID
# %a: Array Task ID
# ----------------------------------------------------------
#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=80G
#SBATCH --time=23:59:00                          # generous; lower once you know per-pixel runtime
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_singlePix_TEST_%A_%a
#SBATCH --output=EMIT_singlePix_TEST_%A_%a.out
#SBATCH --error=EMIT_singlePix_TEST_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1001                             # single task -> one pixel

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
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

cd /projects/anbu8374/


# Load MATLAB
module load matlab/R2024b


# ----------------------------------------------------------
# *** DIRECTORY CONTAINING PER-PIXEL .mat FILES (no trailing slash) ***
# ----------------------------------------------------------
INPUT_DIR="/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver2"

# ----------------------------------------------------------
# *** TEST OUTPUT DIRECTORY (trailing slash required) ***
# Kept separate from take_13/ so the test output does not mix with the
# production run. Change to the real take_13/ dir if you prefer.
# ----------------------------------------------------------
OUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/Droplet_profile_retrievals/take_13_test/"
mkdir -p "${OUT_DIR}"


# ----------------------------------------------------------
# *** WHICH .mat FILE TO TEST ***
# Index into the sorted list of .mat files (0-based). 0 = first file.
# Set to a specific index if you want to test a particular pixel.
# ----------------------------------------------------------
TEST_FILE_INDEX=0


# Get sorted list of .mat files (full paths)
mapfile -t ALL_MAT_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f | sort)
TOTAL_FILES=${#ALL_MAT_FILES[@]}

if [ ${TOTAL_FILES} -eq 0 ]; then
    echo "ERROR: No .mat files found in ${INPUT_DIR}"
    exit 1
fi

if [ ${TEST_FILE_INDEX} -ge ${TOTAL_FILES} ]; then
    echo "ERROR: TEST_FILE_INDEX=${TEST_FILE_INDEX} but only ${TOTAL_FILES} files found."
    exit 1
fi

CURRENT_MAT_FILE=${ALL_MAT_FILES[$TEST_FILE_INDEX]}

# folder_extension_number gives unique INP_OUT/wc/atmmod/mie dirs for this task
FOLDER_EXT_NUM=${SLURM_ARRAY_TASK_ID}


# Give the task its own isolated MATLAB preferences and job directory
TMPDIR=/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "${TMPDIR}"
export TMPDIR


echo " "
echo "=============================================="
echo "EMIT single-pixel SMOKE TEST"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "PRE_MATLAB_LD_LIBRARY_PATH: ${PRE_MATLAB_LD_LIBRARY_PATH}"
echo "TMPDIR: ${TMPDIR}"
echo "Total .mat files found: ${TOTAL_FILES}"
echo "Testing file index ${TEST_FILE_INDEX}: ${CURRENT_MAT_FILE}"
echo "folder_extension_number: ${FOLDER_EXT_NUM}"
echo "Output dir: ${OUT_DIR}"
echo "Starting MATLAB at $(date)"
echo "=============================================="


time matlab -nodesktop -nodisplay -r "\
    addpath(genpath('/projects/anbu8374/Matlab-Research')); \
    addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); \
    addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); \
    clear variables; \
    addLibRadTran_paths; \
    print_status_updates = true; \
    print_libRadtran_err = true; \
    delete_inp_out = true; \
    mat_file_path = '${CURRENT_MAT_FILE}'; \
    folder_extension_number = ${FOLDER_EXT_NUM}; \
    output_dir = '${OUT_DIR}'; \
    [GN_inputs, GN_outputs, tblut_retrieval, acpw_retrieval, folder_paths] = \
        run_retrieval_singlePixel_EMIT_Aqua(mat_file_path, folder_extension_number, \
        print_status_updates, print_libRadtran_err, output_dir, delete_inp_out); \
    exit"

MATLAB_EXIT=$?

echo " "
echo "MATLAB exited with status ${MATLAB_EXIT} at $(date)"
echo "Output .mat files now in ${OUT_DIR}:"
ls -lh "${OUT_DIR}" 2>/dev/null || echo "  (none -- check the .err log)"


# ----------------------------------------------------------
# Cleanup (same as production)
# ----------------------------------------------------------
echo "Cleaning MATLAB temp directory ${TMPDIR}"
rm -rf "${TMPDIR}"

# Remove THIS test's unique libRadtran dirs immediately (single known number)
echo "Removing this test's libRadtran scratch/projects dirs (ext ${FOLDER_EXT_NUM})"
rm -rf "/scratch/alpine/${USER}/EMIT/INP_OUT_${FOLDER_EXT_NUM}"               2>/dev/null || true
rm -rf "/scratch/alpine/${USER}/Mie_Calculations/emit_${FOLDER_EXT_NUM}"      2>/dev/null || true
rm -rf "/projects/${USER}/software/libRadtran-2.0.5/data/wc_${FOLDER_EXT_NUM}"     2>/dev/null || true
rm -rf "/projects/${USER}/software/libRadtran-2.0.5/data/atmmod_${FOLDER_EXT_NUM}" 2>/dev/null || true

echo " "
echo "=============================================="
echo "TEST task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
echo "If MATLAB exit status is 0 and an output .mat exists above, the pipeline works."
echo "=============================================="
