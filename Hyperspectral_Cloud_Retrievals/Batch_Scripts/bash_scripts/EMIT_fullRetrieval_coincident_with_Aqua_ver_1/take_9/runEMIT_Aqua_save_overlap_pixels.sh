#!/bin/bash

# SLURM Script to find overlapping EMIT-Aqua pixels and save per-pixel .mat files
#
# This is the PREPROCESSING step that runs before the retrieval job array.
# It calls save_overlap_data_perPixel_EMIT_Aqua.m for each coincident data
# subdirectory, saving one .mat file per EMIT pixel. Those .mat files are
# then processed by the retrieval job array (runEMIT_Aqua_singlePixel_retrieval.sh).
#
# The coincident_dataPath and sub_directories are hard-coded for the CU Boulder
# Alpine supercomputer. The overlap criteria can be changed in the section below.

# %A: Job ID

#SBATCH --account=ucb762_asc1                   # Ascent Allocation on Alpine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=90G
#SBATCH --time=23:59:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=EMIT_saveOverlap_%A
#SBATCH --output=EMIT_saveOverlap_%A.out
#SBATCH --error=EMIT_saveOverlap_%A.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL

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


# ===========================================================
# *** OVERLAP CRITERIA - EDIT THESE AS NEEDED ***
# ===========================================================
CLD_PHASE="water"
CLD_CVR=1
CLD_TAU_MIN=3
CLD_TAU_MAX=50
H=2
FIND_N_SMALLEST_H=true
H_N_SMALLEST=10
EMIT_PIXELS_PER_MODIS=10
PRINT_STATUS_UPDATES=true
# ===========================================================




# ---------------------------------------------------------------
# *** DEFINE THE DIRECTORY WHERE ALL .mat FILES WILL BE SAVED ***
# *** MUST HAVE TRAILING SLASH '/' AT THE END             ***
# ---------------------------------------------------------------
OUT_DIR="/scratch/alpine/anbu8374/EMIT_pix_overlap_with_Aqua_paper2_ver2/"
# ----------------------------------------------------------


# Start of the job
echo " "
echo "Starting EMIT overlap pixel preprocessing at $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo " "
echo "Criteria:"
echo "  cld_phase = ${CLD_PHASE}"
echo "  cld_cvr = ${CLD_CVR}"
echo "  cld_tau_min = ${CLD_TAU_MIN}"
echo "  cld_tau_max = ${CLD_TAU_MAX}"
echo "  H = ${H}"
echo "  findN_smallest_H = ${FIND_N_SMALLEST_H}"
echo "  H_N_smallest = ${H_N_SMALLEST}"
echo "  emit_pixels_per_modis = ${EMIT_PIXELS_PER_MODIS}"


# Clean MATLAB temp directories
echo " "
echo "Cleaning MATLAB temp directories"
rm -rf ~/.matlab/local_cluster_jobs/R2024b/Job*
rm -rf /tmp/mathworks_${USER}_*


# Run MATLAB
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "\
    addpath(genpath('/projects/anbu8374/Matlab-Research')); \
    clear variables; \
    addLibRadTran_paths; \
    folder_paths = define_EMIT_dataPath_and_saveFolders(5); \
    folder_paths.coincident_dataPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'; \
    sub_directories = {'2023_9_16_T191118_1/', '2023_9_16_T191130_1/', '2023_9_16_T191142_1/', \
                       '2024_1_13_T194658_1/', '2024_5_17_T183906_1/', '2024_5_17_T183918_1/', \
                       '2024_5_17_T183930_1/'}; \
    folder_paths.output_dir = '${OUT_DIR}'; \
    mkdir(folder_paths.output_dir); \
    criteria.cld_phase = '${CLD_PHASE}'; \
    criteria.cld_cvr = ${CLD_CVR}; \
    criteria.cld_tau_min = ${CLD_TAU_MIN}; \
    criteria.cld_tau_max = ${CLD_TAU_MAX}; \
    criteria.H = ${H}; \
    criteria.findN_smallest_H = ${FIND_N_SMALLEST_H}; \
    criteria.H_N_smallest = ${H_N_SMALLEST}; \
    emit_pixels_per_modis = ${EMIT_PIXELS_PER_MODIS}; \
    print_status_updates = ${PRINT_STATUS_UPDATES}; \
    for nn = 1:length(sub_directories), \
        n_saved = save_overlap_data_perPixel_EMIT_Aqua(criteria, folder_paths, print_status_updates, emit_pixels_per_modis, sub_directories{nn}); \
        disp(['Subdirectory ', num2str(nn), ': saved ', num2str(n_saved), ' pixels to ', folder_paths.output_dir]); \
    end; \
    exit"


echo " "
echo "Finished EMIT overlap pixel preprocessing at $(date)"
