#!/bin/bash

# Resubmit the 12 jobs that exceeded the 24h time limit on amilan normal.
# These are submitted to the amilan 'long' partition (qos=long), which allows
# up to 7 days of wall time and a maximum of 150 concurrent jobs per user.
#
# Failed tasks and their measurement indices:
#   sza23 (offset=200): task 239 -> meas 39
#   sza33 (offset=300): task 335 -> meas 35
#   sza47 (offset=500): tasks 526,533,544,550,553 -> meas 26,33,44,50,53
#   sza54 (offset=600): tasks 628,632,641,647,648 -> meas 28,32,41,47,48
#
# Submit each block independently with:
#   sbatch rerun_failed_long_partition.sh <sza_block>
# where <sza_block> is one of: sza23, sza33, sza47, sza54
#
# Or submit all four blocks at once:
#   for block in sza23 sza33 sza47 sza54; do sbatch rerun_failed_long_partition.sh $block; done

# -------------------------------------------------------
# Read which SZA block to run from the first argument
# -------------------------------------------------------
SZA_BLOCK=${1:-sza23}

# -------------------------------------------------------
# Configure per-block settings
# -------------------------------------------------------
case "$SZA_BLOCK" in
    sza23)
        SZA=23.4343
        offset=200
        ARRAY_IDS="239"
        JOB_NAME="rerun_long_sza23"
        ;;
    sza33)
        SZA=33.3806
        offset=300
        ARRAY_IDS="335"
        JOB_NAME="rerun_long_sza33"
        ;;
    sza47)
        SZA=47.9277
        offset=500
        ARRAY_IDS="526,533,544,550,553"
        JOB_NAME="rerun_long_sza47"
        ;;
    sza54)
        SZA=54.0142
        offset=600
        ARRAY_IDS="628,632,641,647,648"
        JOB_NAME="rerun_long_sza54"
        ;;
    *)
        echo "ERROR: Unknown SZA block '$SZA_BLOCK'. Choose one of: sza23, sza33, sza47, sza54"
        exit 1
        ;;
esac

# -------------------------------------------------------
# Resubmit via sbatch with all settings inlined
# -------------------------------------------------------
sbatch <<EOF
#!/bin/bash
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --mem=90G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${JOB_NAME}_%A_%a.out
#SBATCH --error=${JOB_NAME}_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=${ARRAY_IDS}

# Load modules
ml purge
ml gcc/11.2.0
ml netcdf/4.8.1
ml perl/5.36.0
ml texlive/2021

# Set up environment paths
export PATH=/projects/\$USER/software/libRadtran-2.0.5/:$PATH
export PATH=/projects/\$USER/software/libRadtran-2.0.5/data/:$PATH
export PATH=/projects/\$USER/software/libRadtran-2.0.5/bin/:$PATH
export GSL_BIN=/projects/\$USER/software/gsl-2.6/bin
export GSL_LIB=/projects/\$USER/software/gsl-2.6/lib
export GSL_INC=/projects/\$USER/software/gsl-2.6/include
export LD_LIBRARY_PATH=\$GSL_LIB:\$LD_LIBRARY_PATH
export INSTALL_DIR=/projects/\$USER/software/libRadtran-2.0.5
export PATH=\$GSL_BIN:\$PATH

# Capture the correct LD_LIBRARY_PATH before MATLAB contaminates it
export PRE_MATLAB_LD_LIBRARY_PATH=\$LD_LIBRARY_PATH

cd /projects/anbu8374/

module load matlab/R2024b

SZA=${SZA}
offset=${offset}
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/dataSet_created_on_31_March_2026/"
mkdir -p \$output_dir

export TMPDIR=/scratch/alpine/\${USER}/matlab_tmp_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}
mkdir -p \$TMPDIR

# Stagger: only a handful of tasks per block so a short stagger is sufficient
sleep \$((\${SLURM_ARRAY_TASK_ID} % 60))

echo "Starting MATLAB job for measurement \${SLURM_ARRAY_TASK_ID} at \$(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(\${SLURM_ARRAY_TASK_ID}); measurement_idx = \${SLURM_ARRAY_TASK_ID} - \${offset}; hysics_refl_from_vocals_and_era5_SZA_loopGeometry_ver3(folder_paths, measurement_idx, \${SZA}, '\${output_dir}'); exit"

echo "Finished MATLAB job for measurement \${SLURM_ARRAY_TASK_ID} at \$(date)"

echo "Removing TMPDIR: \$TMPDIR"
rm -rf "\$TMPDIR"

echo "Removing INP_OUT directory for task \${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/\${USER}/HySICS/INP_OUT_\${SLURM_ARRAY_TASK_ID}"

echo "Removing atmmod and wc directories for task \${SLURM_ARRAY_TASK_ID}..."
rm -rf "/scratch/alpine/\${USER}/software/libRadtran-2.0.5/data/atmmod_\${SLURM_ARRAY_TASK_ID}"
rm -rf "/scratch/alpine/\${USER}/software/libRadtran-2.0.5/data/wc_\${SLURM_ARRAY_TASK_ID}"

echo "Pruning scratch directories older than 7 days..."
find /scratch/alpine/\${USER}/                  -maxdepth 1 -name "matlab_tmp_*"       -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/\${USER}/HySICS/           -maxdepth 1 -name "INP_OUT_*"          -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
find /scratch/alpine/\${USER}/Mie_Calculations/ -maxdepth 1 -name "Mie_Calculations_*" -type d -mtime +7 -exec rm -rf {} \; 2>/dev/null || true
EOF
