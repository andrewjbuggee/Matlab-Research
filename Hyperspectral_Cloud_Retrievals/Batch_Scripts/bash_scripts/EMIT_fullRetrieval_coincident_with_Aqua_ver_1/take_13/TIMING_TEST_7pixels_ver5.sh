#!/bin/bash

# ======================================================================
# TIMING TEST — 7 EMIT-Aqua single-pixel retrievals, ONE PER ARRAY TASK,
# all running concurrently, to measure per-pixel wall time with the current
# settings (225 bands, 16 DISORT streams, array_length_newMax=15, forward-model
# Jacobian computed once at the initial guess).
#
# Submit AFTER take_9 preprocessing has written the 7 per-pixel .mat files:
#   sbatch TIMING_TEST_7pixels_ver5.sh
#
# Each task processes exactly one pixel. Because all 7 run at once (no %throttle),
# total turnaround is roughly the slowest single pixel. Read the per-pixel wall
# time from the `real ...` line in each task's .err (or the Start/Finish stamps
# in each .out) -- see the summary one-liner at the bottom of this file.
#
# Resources match the intended production run (40 cpus / 100G) so the timing is
# representative.
# ======================================================================

# %A: Job ID   %a: Array Task ID
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --time=23:59:00                          # normal-QOS max; a pixel >24h would be killed (and that itself is useful info)
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_timingTest_%A_%a
#SBATCH --output=EMIT_timingTest_%A_%a.out
#SBATCH --error=EMIT_timingTest_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-7                              # one task per pixel, all concurrent. Set 1-N if take_9 saved N != 7 files.

# Load modules
ml purge
ml gcc/11.2.0
ml netcdf/4.8.1
ml perl/5.36.0
ml texlive/2021

# libRadtran + GSL
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

cd /projects/anbu8374/

# Load MATLAB
module load matlab/R2024b


# ----------------------------------------------------------
# Per-pixel .mat files written by take_9 (no trailing slash)
# ----------------------------------------------------------
INPUT_DIR="/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver2"

# ----------------------------------------------------------
# Test output directory (persistent; kept separate from production take_13)
# ----------------------------------------------------------
OUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/Droplet_profile_retrievals/take_13_timingTest/"
mkdir -p "${OUT_DIR}"


# ----------------------------------------------------------
# Per-task node-local scratch for parpool storage + libRadtran I/O, with a
# Lustre fallback; and a per-task MATLAB preferences dir.
# ----------------------------------------------------------
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
export MATLAB_PREFDIR="${TMPDIR}/matlab_prefs/R2024b"
mkdir -p "${MATLAB_PREFDIR}"
echo "MATLAB_PREFDIR: ${MATLAB_PREFDIR}"


# ----------------------------------------------------------
# Pick this task's pixel: the (task_id)-th .mat file (sorted)
# ----------------------------------------------------------
mapfile -t ALL_MAT_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f | sort)
TOTAL_FILES=${#ALL_MAT_FILES[@]}
FILE_IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))          # 0-based index

if [ ${TOTAL_FILES} -eq 0 ]; then
    echo "ERROR: no .mat files in ${INPUT_DIR} -- run take_9 preprocessing first."
    exit 1
fi
if [ ${FILE_IDX} -ge ${TOTAL_FILES} ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID}: only ${TOTAL_FILES} files found; nothing to do."
    exit 0
fi

CURRENT_MAT_FILE=${ALL_MAT_FILES[$FILE_IDX]}
FOLDER_EXT_NUM=${SLURM_ARRAY_TASK_ID}

echo " "
echo "=============================================="
echo "EMIT TIMING TEST  task ${SLURM_ARRAY_TASK_ID} of ${TOTAL_FILES}"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "Pixel file: ${CURRENT_MAT_FILE}"
echo "Starting MATLAB at $(date)"
echo "=============================================="


# The `time` wrapper writes the per-pixel wall time (real/user/sys) to the .err.
time matlab -nodesktop -nodisplay -r "\
    try; \
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
        exit(0); \
    catch ME; \
        fprintf(2, '\n[RETRIEVAL FAILED] %s: %s\n', ME.identifier, ME.message); \
        for k = 1:numel(ME.stack); fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; \
        exit(1); \
    end"
MATLAB_STATUS=$?

echo " "
echo "Task ${SLURM_ARRAY_TASK_ID} finished: MATLAB status ${MATLAB_STATUS} at $(date)"


# ----------------------------------------------------------
# Cleanup (INP_OUT/wc/atmmod live under $TMPDIR via define_EMIT; Mie on /scratch)
# ----------------------------------------------------------
rm -rf "${TMPDIR}"
rm -rf "/scratch/alpine/${USER}/EMIT/INP_OUT_${FOLDER_EXT_NUM}"          2>/dev/null || true
rm -rf "/scratch/alpine/${USER}/Mie_Calculations/emit_${FOLDER_EXT_NUM}" 2>/dev/null || true

exit ${MATLAB_STATUS}

# ----------------------------------------------------------
# After all 7 finish, summarize the per-pixel wall times:
#   grep -H "real" EMIT_timingTest_*.err
# or as a quick table (minutes):
#   for f in EMIT_timingTest_*.err; do \
#     t=$(grep -oE 'real[[:space:]]+[0-9]+m[0-9.]+s' "$f"); echo "$f  $t"; done
# ----------------------------------------------------------
