#!/bin/bash
INPUT_DIR="/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra"
mapfile -t ALL_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.mat" -type f -printf "%f\n" | sort)

TOTAL_FILES=${#ALL_FILES[@]}
ARRAY_SIZE=36
FILES_PER_JOB=$(( (TOTAL_FILES + ARRAY_SIZE - 1) / ARRAY_SIZE ))

echo "Total files: ${TOTAL_FILES}"
echo "Array jobs: ${ARRAY_SIZE}"  
echo "Files per job: ${FILES_PER_JOB}"
echo "Total CPU cores needed: $((ARRAY_SIZE * 16))"
echo "Nodes needed: $(( (ARRAY_SIZE * 16 + 63) / 64 ))"