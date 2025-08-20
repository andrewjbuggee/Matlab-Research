#!/bin/bash

#

#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=amilan
#SBATCH --qos=mem
#SBATCH --mem=215G
#SBATCH --ntasks=10
#SBATCH --job-name=test_retrieval_HySICS_rTop_10_rBot_tauC_11_tcwv_14_vza_7_test_highMem
#SBATCH --output=test_retrieval_HySICS_rTop_10_rBot_tauC_11_tcwv_14_vza_7_highMem.out
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

module load matlab

echo "Starting MATLAB job at $(date)"

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); test_retrieval_HySICS_rTop_10_rBot_tauC_11_tcwv_14; exit"

echo "Finished MATLAB job at $(date)"
