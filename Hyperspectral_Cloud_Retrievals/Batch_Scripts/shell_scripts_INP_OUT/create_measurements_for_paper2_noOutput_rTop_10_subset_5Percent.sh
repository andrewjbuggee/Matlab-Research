#!/bin/bash

#

#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=17G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=create_measurements_for_paper2_noOutput_rTop_10_subset_5Percent
#SBATCH --output=measurements_for_paper2_rTop_10_subset_5Percent.out
#SBATCH --error=create_measurements_for_paper2_noOutput_rTop_10_subset_5Percent.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL


ml purge
ml gcc/11.2.0
ml netcdf/4.8.1
ml perl/5.36.0
ml texlive/2021


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


# Define the directory where the files generated should be saved
# ----------------------------------------------------------
# *** MODIFY THIS DIRECTORY BASED ON THE LOCATION OF THE MEASUREMENTS ***
# *** CANNOT HAVE TRAILING SLASH '/' AT THE END         ***
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/paper2_variableSweep/rTop_10/vza_7_vaz_210_sza_10_saz_91_subset/"



module load matlab

# Run MATLAB once with all files for this job
echo " "
echo "Starting MATLAB at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; folder_paths = define_folderPaths_for_HySICS('${SLURM_ARRAY_JOB_ID}'); folder_paths.HySICS_simulated_spectra = '${INPUT_DIR}/'; print_libRadtran_err = true; generate_hysics_measurements_noOutput_paper2_rTop_10_5Percent; exit"

echo " "
echo "Finished MATLAB job array task ${SLURM_ARRAY_TASK_ID} at $(date)"
