#!/bin/bash

# Hybrid approach: Process multiple files per job to maximize node utilization
# In the directory below, there are 12 files
# With 12 jobs and 12 files, each job will process 1 file

# SLURM Job Array Script to run MATLAB retrievals on multiple files in parallel
# This will run the same analysis on multiple files within a specified directory

# Because this is a job array, each job will request resources independently
# This means each job will request N ntasks, on N nodes with N cpus-per-task

# %A: Job ID
# %a: Array Task ID

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=36G
#SBATCH --time=3:00:00     # Longer time for multiple files
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=full_retrieval_hysics_rTop_10_subset_pt1Percent_%A_%a
#SBATCH --output=full_retrieval_hysics_rTop_10_subset_pt1Percent_%A_%a.out
#SBATCH --error=full_retrieval_hysics_rTop_10_subset_pt1Percent_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-12       # 12 jobs × 1 file each = 12 files

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
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset"

RETRIEVED_PROFS_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/"

# Get list of all files that have 0.1% measurement uncertainty
mapfile -t ALL_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "simulated_spectra_HySICS_reflectance_66bands_0.1%*.mat" -type f -printf "%f\n" | sort)

# Calculate which files this job should process
FILES_PER_JOB=1
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


# -------------------------------------------------------------
# Debugging section to check is file array is defined correctly
echo " "
echo "=== DEBUG INFO ==="
echo "SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Total files found: ${#ALL_FILES[@]}"
echo "START_IDX: ${START_IDX}"
echo "END_IDX: ${END_IDX}"
echo "FILE_ARRAY contents: [${FILE_ARRAY}]"

# List the first few files found:
echo "First few files found:"
for (( j=0; j<5 && j<${#ALL_FILES[@]}; j++ )); do
    echo "  [$j]: ${ALL_FILES[$j]}"
done
echo "=================="

# Check if FILE_ARRAY is empty before running MATLAB
if [ -z "$FILE_ARRAY" ]; then
    echo "ERROR: No files assigned to this job!"
    echo "This usually means the array indexing is wrong."
    exit 1
fi
# ----------------------------------------------

# Run MATLAB once with all files for this job
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS('${SLURM_ARRAY_TASK_ID}'); folder_paths.HySICS_simulated_spectra = '${INPUT_DIR}/'; folder_paths.HySICS_retrievals = '${RETRIEVED_PROFS_DIR}'; print_status_updates = false; print_libRadtran_err = false; file_list = {${FILE_ARRAY}}; [tblut_retrieval, acpw_retrieval, GN_inputs, GN_outputs] = run_retrieval_dropletProfile_HySICS_ver3(file_list, folder_paths, print_status_updates, print_libRadtran_err); exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"