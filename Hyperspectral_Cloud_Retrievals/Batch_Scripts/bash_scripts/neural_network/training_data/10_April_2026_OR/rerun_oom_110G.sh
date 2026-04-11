#!/bin/bash

# Resubmit the 4 ORACLES measurement indices that were OOM-killed at 90G.
# Measurement indices 132, 154, 206, and 234 are consistently OOM-killed
# across all SZA scripts. This script resubmits them at 110G.
#
# These measurement indices correspond to the following task IDs per SZA:
#   sza0  (offset=100):  tasks 232, 254, 306, 334
#   sza23 (offset=400):  tasks 532, 554, 606, 634
#   sza33 (offset=700):  tasks 832, 854, 906, 934
#   sza41 (offset=1000): tasks 1132, 1154, 1206, 1234
#   sza47 (offset=1300): tasks 1432, 1454, 1506, 1534
#   sza54 (offset=1600): tasks 1732, 1754, 1806, 1834
#   sza59 (offset=2000): tasks 2132, 2154, 2206, 2234
#   sza65 (offset=2300): tasks 2432, 2454, 2506, 2534
#
# Usage: sbatch rerun_oom_110G.sh <sza_block>
# where <sza_block> is one of: sza0, sza23, sza33, sza41, sza47, sza54, sza59, sza65
#
# Submit all blocks at once:
#   for block in sza0 sza23 sza33 sza41 sza47 sza54 sza59 sza65; do
#       sbatch rerun_oom_110G.sh $block
#   done

SZA_BLOCK=${1:-sza0}

case "$SZA_BLOCK" in
    sza0)
        SZA=0
        offset=100
        ARRAY_IDS="232,254,306,334"
        ;;
    sza23)
        SZA=23.4343
        offset=400
        ARRAY_IDS="532,554,606,634"
        ;;
    sza33)
        SZA=33.3806
        offset=700
        ARRAY_IDS="832,854,906,934"
        ;;
    sza41)
        SZA=41.1882
        offset=1000
        ARRAY_IDS="1132,1154,1206,1234"
        ;;
    sza47)
        SZA=47.9277
        offset=1300
        ARRAY_IDS="1432,1454,1506,1534"
        ;;
    sza54)
        SZA=54.0142
        offset=1600
        ARRAY_IDS="1732,1754,1806,1834"
        ;;
    sza59)
        SZA=59.6619
        offset=2000
        ARRAY_IDS="2132,2154,2206,2234"
        ;;
    sza65)
        SZA=65.0000
        offset=2300
        ARRAY_IDS="2432,2454,2506,2534"
        ;;
    *)
        echo "ERROR: Unknown SZA block '$SZA_BLOCK'."
        echo "Choose one of: sza0, sza23, sza33, sza41, sza47, sza54, sza59, sza65"
        exit 1
        ;;
esac

sbatch <<EOF
#!/bin/bash
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=23:59:59
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=110G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=rerun_oom_OR_${SZA_BLOCK}
#SBATCH --output=rerun_oom_OR_${SZA_BLOCK}_%A_%a.out
#SBATCH --error=rerun_oom_OR_${SZA_BLOCK}_%A_%a.err
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

export PRE_MATLAB_LD_LIBRARY_PATH=\$LD_LIBRARY_PATH

cd /projects/anbu8374/

module load matlab/R2024b

SZA=${SZA}
offset=${offset}
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/dataSet_created_on_10_April_2026_OR/"
mkdir -p \$output_dir

export TMPDIR=/scratch/alpine/\${USER}/matlab_tmp_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}
mkdir -p \$TMPDIR

# Only 4 tasks per block, minimal stagger needed
sleep \$((\${SLURM_ARRAY_TASK_ID} % 60))

echo "Starting MATLAB job for measurement \${SLURM_ARRAY_TASK_ID} at \$(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS(\${SLURM_ARRAY_TASK_ID}); measurement_idx = \${SLURM_ARRAY_TASK_ID} - \${offset}; hysics_refl_from_oracles_and_era5_SZA_loopGeometry(folder_paths, measurement_idx, \${SZA}, '\${output_dir}'); exit"

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
