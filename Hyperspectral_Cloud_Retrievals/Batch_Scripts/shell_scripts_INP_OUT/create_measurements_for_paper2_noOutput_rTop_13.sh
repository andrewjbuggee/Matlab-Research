#!/bin/bash

#

#SBATCH --nodes=1
#SBATCH --time=15:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=17G
#SBATCH --ntasks=10
#SBATCH --job-name=create_measurements_for_paper2_noOutput_rTop_13
#SBATCH --output=measurements_for_paper2_rTop_13_newGeo_rev4.out
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

time matlab -nodesktop -nodisplay -r "addpath(genpath('/projects/anbu8374/Matlab-Research')); addpath(genpath('/scratch/alpine/anbu8374/HySICS/INP_OUT/')); addpath(genpath('/scratch/alpine/anbu8374/Mie_Calculations/')); generate_hysics_measurements_noOutput_paper2_rTop_13; exit"

echo "Finished MATLAB job at $(date)"
