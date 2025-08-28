#!/bin/bash

# Testing job_array


#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --array=1-2
#SBATCH --mem=35G
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test_retrieval_array_amilan_10Tasks_%A_%a
#SBATCH --output=test_retrieval_array_amilan_10Tasks_%A_%a.out
#SBATCH --error=test_retrieval_array_amilan_10Tasks_%A_%a.err
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

module load matlab/R2024b

echo "Starting MATLAB job at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); clear variables; addLibRadTran_paths; [folder_paths, which_computer] = define_folderPaths_for_HySICS(3); print_status_updates = true; print_libRadtran_err = true; test_retrieval_HySICS_no_msgs_orDefinedFolder; exit"

echo "Finished MATLAB job at $(date)"
