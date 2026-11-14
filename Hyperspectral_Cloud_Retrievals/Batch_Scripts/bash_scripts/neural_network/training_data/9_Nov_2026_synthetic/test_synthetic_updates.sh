#!/bin/bash
#
# TEST run of the ver2 synthetic-cloud RT pipeline (14 Nov 2026) on Alpine's
# testing QoS: run_synthetic_spectra_chunk.m -> hysics_refl_from_synthetic_NN_inputs_ver2.m
# (4 uvspec runs per cloud instead of 636 per-band runs; parfor over clouds;
# one output .mat for the whole task).
#
# Clouds: the five checked earlier, then clouds 1-19, so 24 clouds fill the 12
# workers twice and the log reports a median time per cloud on Alpine
# ("=== CHUNK SUMMARY ... median N s per cloud"), which sets CHUNK_SIZE in
# synthetic_12_Nov_2026.sh:
#     553    typical cloud      (tau_c = 7.8, the median)
#   15022    edge case          LWC = 0 at the cloud-top level
#   27045    edge case          cloud-top LWC = 8.97 g/m^3 (largest in the set)
#    1105    thinnest cloud     tau_c = 3.00
#   21014    thickest cloud     tau_c = 39.97
# Override:  CLOUD_IDS_OVERRIDE="553 15022" sbatch test_synthetic_updates.sh
#
# BEFORE SUBMITTING (on Alpine):
#   1. git pull in /projects/anbu8374/Matlab-Research (needs
#      run_synthetic_spectra_chunk.m; the pre-flight block checks)
#   2. mkdir -p /projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch_test
#   3. module load slurm/alpine, then  sbatch test_synthetic_updates.sh
#
# Resources: cpu-testing allows 1 h and up to 16 cores (CURC alpine-hardware
# page). Each worker needs ~4 GB (the largest uvspec segment run peaked at
# 3.1 GB on the laptop, plus the MATLAB worker), so 12 workers + client is
# ~52 GB; --mem=56G is billed as ~15 cores, inside the 16-core limit.
# The 13-Nov run lost 77 min to addpath(genpath('/scratch/.../Mie_Calculations/'))
# (2,885 old working folders); that call is gone, so start-up should take
# minutes. Staging the Mie tables took 561 s from /scratch on 13 Nov.
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=acpu        # Alpine CPU partition (Nov 2026 sinfo lists no atesting or amilan)
#SBATCH --qos=cpu-testing      # 1 h, up to 16 cores
#SBATCH --mem=56G               # 12 workers x ~4 GB + client; billed as ~15 cores (cpu-testing max 16)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=synthetic_spectra_ver2_TEST
#SBATCH --chdir=/projects/anbu8374/slurm_logs/12_Nov_2026_synthetic_batch_test
#SBATCH --output=create_meas_synthetic_NN_trainingData_TEST_12Nov_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_TEST_12Nov_%A_%a.err
#SBATCH --mail-user=anbu8374@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1

# --- Test parameters -----------------------------------------------------
CLOUD_IDS=${CLOUD_IDS_OVERRIDE:-"553 15022 27045 1105 21014 $(seq -s ' ' 1 19)"}
N_SEGMENTS=${N_SEGMENTS_OVERRIDE:-4}     # uvspec runs per cloud
FILE_TAG="_TEST"
N_TOTAL=${N_TOTAL_OVERRIDE:-100000}      # length of the 'cloud' dimension in the input .nc
# -------------------------------------------------------------------------

# Modules
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
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

cd /projects/anbu8374/ || exit 1

module load matlab/R2024b

# Input file (the .nc produced by 07_build_training_inputs.py on 12 Nov 2026)
input_file="/scratch/alpine/anbu8374/neural_network_training_data/training_inputs_jointMVN_N100000_L7_2026-11-12_ver2.nc"

# Test outputs go to their own folder so they never mix with a production set.
# (Trailing slash required by the MATLAB driver.)
output_dir="/scratch/alpine/anbu8374/neural_network_training_data/synthetic_dataSet_TEST_12_Nov_2026/"
mkdir -p "$output_dir"

# --- Pre-flight checks: fail in seconds, not after MATLAB starts -----------
code_dir="/projects/anbu8374/Matlab-Research"
echo "Matlab-Research commit: $(git -C "$code_dir" log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo 'unknown (not a git checkout?)')"
if [ -z "$(find "$code_dir" -name write_wc_file_mixed_layers.m -print -quit)" ]; then
    echo "ERROR: write_wc_file_mixed_layers.m not found under $code_dir. git pull first." >&2
    exit 1
fi
if [ ! -f "$input_file" ]; then
    echo "ERROR: input file not found: $input_file (scp it from the laptop first)." >&2
    exit 1
fi
if [ -z "$(find "$code_dir" -name run_synthetic_spectra_chunk.m -print -quit)" ]; then
    echo "ERROR: run_synthetic_spectra_chunk.m not found under $code_dir. git pull first." >&2
    exit 1
fi
# Guard against the pre-azimuth-fix input file: the corrected
# training_inputs_jointMVN_N100000_L7_2026-11-12_ver2.nc (07, 12 Nov 2026,
# 22:04) has this MD5; the old wrong-azimuth _ver2 was d8e67d5b5db5157e9ff565fc958d0f1e.
EXPECTED_MD5=${EXPECTED_MD5_OVERRIDE:-f619a48cffea9631ea3899faa1166af4}
actual_md5=$(md5sum "$input_file" | cut -d' ' -f1)
if [ "$actual_md5" != "$EXPECTED_MD5" ]; then
    echo "ERROR: $input_file has MD5 $actual_md5, expected $EXPECTED_MD5 (azimuth-corrected file)." >&2
    exit 1
fi
echo "Input MD5 OK: $actual_md5"
for id in $CLOUD_IDS; do
    if [ "$id" -lt 1 ] || [ "$id" -gt "$N_TOTAL" ]; then
        echo "ERROR: cloud id $id is outside 1..$N_TOTAL" >&2
        exit 1
    fi
done
# MATLAB vector literal, e.g. "553 15022 27045" -> "[553 15022 27045]"
cloud_id_vec="[${CLOUD_IDS}]"
echo "Task ${SLURM_ARRAY_TASK_ID}: test clouds ${cloud_id_vec}"
echo "Input : $input_file"
echo "Output: $output_dir"
# -------------------------------------------------------------------------

# Per-task scratch on node-local /tmp instead of Lustre (keeps parpool's
# JobStorageLocation off the shared filesystem). Falls back to Lustre if /tmp
# is unwritable.
# CURC recommends $SLURM_SCRATCH: the job's own directory on the node's local
# SSD (~300 GB/node), deleted automatically when the job ends
# (curc.readthedocs.io, compute/filesystems.html, read 13 Nov 2026).
# Falls back to /tmp if it is not set.
TMPDIR_LOCAL="${SLURM_SCRATCH:-/tmp}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if mkdir -p "$TMPDIR_LOCAL" 2>/dev/null && [ -w "$TMPDIR_LOCAL" ]; then
    export TMPDIR="$TMPDIR_LOCAL"
    echo "TMPDIR (node-local, SLURM_SCRATCH=${SLURM_SCRATCH:-unset}): $TMPDIR"
else
    export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
    echo "TMPDIR (Lustre fallback): $TMPDIR"
fi
echo "TMPDIR headroom:"
df -h "$TMPDIR"

# --- Stage the custom Mie tables on node-local disk (added 13 Nov 2026) ----
# Every uvspec run opens the cloud's custom Mie table (~130 MB each, 1-50 um
# family): 4 runs per cloud with the ver2 driver, 636 with ver1. CURC asks
# that /projects not carry intensive job I/O. Copy the 48 tables once per task to $SLURM_SCRATCH (local SSD, deleted
# at job end) and point MATLAB at them via CUSTOM_MIE_TABLES_DIR, which
# find_custom_mieTable_closest_to_alpha_profile.m honours. Source: a /scratch
# copy if one exists (Lustre is built for job I/O), otherwise /projects.
# If the copy is incomplete, fall back to the original tables with a warning.
MIE_PATTERN="wc_mieTable_gamma_rEff_1-50microns_gammaDist_alpha_*.cdf"
MIE_N_EXPECTED=48
MIE_SCRATCH_SRC="/scratch/alpine/${USER}/custom_mieTables"
MIE_PROJECTS_SRC="/projects/${USER}/software/libRadtran-2.0.5/data/wc/custom_mieTables"
MIE_SRC=""
for cand in "$MIE_SCRATCH_SRC" "$MIE_PROJECTS_SRC"; do
    if [ "$(find "$cand" -maxdepth 1 -name "$MIE_PATTERN" 2>/dev/null | wc -l)" -ge "$MIE_N_EXPECTED" ]; then
        MIE_SRC="$cand"
        break
    fi
done
# Say loudly when the /scratch copy is gone (e.g. removed by the 90-day
# /scratch purge) so it can be re-copied; the job still runs from /projects.
if [ "$MIE_SRC" != "$MIE_SCRATCH_SRC" ]; then
    n_scratch=$(find "$MIE_SCRATCH_SRC" -maxdepth 1 -name "$MIE_PATTERN" 2>/dev/null | wc -l | tr -d ' ')
    msg="NOTICE: /scratch Mie-table copy is missing or incomplete (${n_scratch} of ${MIE_N_EXPECTED} files in ${MIE_SCRATCH_SRC}); staging from ${MIE_SRC:-nowhere} instead. Re-copy on a login node with:
    mkdir -p ${MIE_SCRATCH_SRC} && cp ${MIE_PROJECTS_SRC}/${MIE_PATTERN} ${MIE_SCRATCH_SRC}/"
    echo "$msg"          # .out log
    echo "$msg" >&2      # .err log, where it stands out
fi
MIE_LOCAL="${SLURM_SCRATCH:-$TMPDIR}/custom_mieTables"
if [ -n "$MIE_SRC" ] && mkdir -p "$MIE_LOCAL"; then
    t_copy=$(date +%s)
    echo "Mie tables: copying ${MIE_N_EXPECTED} files from $MIE_SRC to $MIE_LOCAL (started $(date '+%H:%M:%S')) ..."
    cp "$MIE_SRC"/$MIE_PATTERN "$MIE_LOCAL"/ 2>/dev/null
    n_copied=$(find "$MIE_LOCAL" -maxdepth 1 -name "$MIE_PATTERN" | wc -l)
    if [ "$n_copied" -ge "$MIE_N_EXPECTED" ]; then
        export CUSTOM_MIE_TABLES_DIR="$MIE_LOCAL/"
        echo "Mie tables: staged $n_copied files ($(du -sh "$MIE_LOCAL" | cut -f1)) from $MIE_SRC to $MIE_LOCAL in $(( $(date +%s) - t_copy )) s"
    else
        echo "WARNING: only $n_copied of $MIE_N_EXPECTED Mie tables copied; libRadtran will read them from /projects" >&2
        rm -rf "$MIE_LOCAL"
    fi
else
    echo "WARNING: no complete Mie-table source found; libRadtran will read the default /projects tables" >&2
fi
# -------------------------------------------------------------------------

# (No start-up stagger: a single test task has nothing to collide with.)

echo "Starting MATLAB test job at $(date)"

# MATLAB preferences in this task's scratch to avoid collisions with /home.
export MATLAB_PREFDIR="$TMPDIR/matlab_prefs/R2024b"
mkdir -p "$MATLAB_PREFDIR"
echo "MATLAB_PREFDIR: $MATLAB_PREFDIR"

# Same MATLAB call as the production array, except: the explicit cloud list,
# and a 60 s parpool retry pause (production waits 300-420 s).
# Start-up steps are timestamped ([time] elapsed-s step): the 13-Nov test spent
# ~53 min between MATLAB launch and the first cloud, so each step is timed.
matlab_status=0
time matlab -nodesktop -nodisplay -r "try, t_su = tic; fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'MATLAB up'); p = strsplit(genpath('/projects/anbu8374/Matlab-Research'), pathsep); p = p(~contains(p, '/.') & ~cellfun(@isempty, p)); addpath(strjoin(p, pathsep)); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'path: Matlab-Research (no hidden folders)'); fprintf('MATLAB path entries: %d\n', numel(strsplit(path, pathsep))); addLibRadTran_paths; fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'addLibRadTran_paths'); folder_paths = define_folderPaths_for_HySICS(${SLURM_ARRAY_TASK_ID}); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'define_folderPaths_for_HySICS'); for pool_attempt = 1:2, try, fprintf('\n=== PARPOOL START ATTEMPT %d of 2 ===\n', pool_attempt); start_parallel_pool(folder_paths.which_computer); break; catch ME_pool, fprintf(2, '\n[PARPOOL START FAILED attempt %d] %s: %s\n', pool_attempt, ME_pool.identifier, ME_pool.message); old_pool = gcp('nocreate'); if ~isempty(old_pool), delete(old_pool); end; if pool_attempt < 2, pause(60); else, rethrow(ME_pool); end; end; end; fprintf('\nWorkers in pool: %d\n', gcp('nocreate').NumWorkers); fprintf('[%s] %7.1f s  %s\n', char(datetime('now')), toc(t_su), 'parpool ready'); chunk_settings = struct('n_segments', ${N_SEGMENTS}, 'file_tag', '${FILE_TAG}'); summ = run_synthetic_spectra_chunk('${input_file}', ${cloud_id_vec}, folder_paths, '${output_dir}', chunk_settings); if summ.n_failed > 0, exit(2); else, exit(0); end; catch ME, fprintf(2, '\n[FATAL MATLAB ERROR] %s: %s\n', ME.identifier, ME.message); for k = 1:numel(ME.stack), fprintf(2, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line); end; exit(1); end"
matlab_status=$?

echo "Finished MATLAB test job at $(date)"

# What the test produced (newest first), so the .out log shows it directly.
echo "Files written to $output_dir:"
ls -lt "$output_dir" | head -20

# Cleanup — INP_OUT, wc and atmmod live under $TMPDIR.
echo "Removing TMPDIR: $TMPDIR"
rm -rf "$TMPDIR"

# (The production script's 7-day scratch pruning is left out of the test.)

if [ "$matlab_status" -ne 0 ]; then
    echo "MATLAB exited with status ${matlab_status}; marking this Slurm task as failed."
    exit "$matlab_status"
fi
