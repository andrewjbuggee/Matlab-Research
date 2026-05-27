# Running libRadtran on CURC Alpine — Working Notes

A consolidated reference covering how `uvspec` is invoked on Alpine, the
SLURM submission pattern that works, MATLAB-side parallelization, the
quirks that have bitten us, and the bottlenecks we've identified and
fixed. Synthesized from troubleshooting sessions during late
April–early May 2026 (synthetic NN training-data generation, 300,001
clouds × 7-level droplet profiles).

Cluster: CURC Alpine. Account `ucb762_asc1`, user `anbu8374`.

---

## 1. libRadtran install and invocation

libRadtran is **user-built under `/projects`**, not a CURC module:

```
/projects/$USER/software/libRadtran-2.0.5/        # INSTALL_DIR
/projects/$USER/software/libRadtran-2.0.5/bin/    # uvspec, mie
/projects/$USER/software/libRadtran-2.0.5/data/   # default data tree
/projects/$USER/software/gsl-2.6/                 # GSL dep, also user-built
```

Every batch script begins with this environment block (verbatim from
`synthetic_4_May_2026.sh`):

```bash
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

# Critical: MATLAB clobbers LD_LIBRARY_PATH on startup. CURC's MATLAB
# wrapper restores from PRE_MATLAB_LD_LIBRARY_PATH if it's set.
export PRE_MATLAB_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

module load matlab/R2024b
```

**File layout per running task** (driven by `define_folderPaths_for_HySICS.m`,
`curc` branch):

| Purpose       | Path                                                                                 | Lives on    |
|---------------|--------------------------------------------------------------------------------------|-------------|
| `.INP`/`.OUT` | `$TMPDIR/HySICS/INP_OUT_<task_id>/`                                                  | node-local `/tmp` |
| `wc` / `atmmod` per-task copies | `$TMPDIR/.../data/wc_<task_id>/`, `.../atmmod_<task_id>/`          | node-local `/tmp` |
| Mie tables    | `/scratch/alpine/$USER/Mie_Calculations/Mie_Calculations_<task_id>/`                 | Lustre — kept on `/scratch` so they can be reused across runs |
| Outputs       | `/scratch/alpine/$USER/neural_network_training_data/<dataset>/`                      | Lustre |

The task-id suffix is essential — concurrent array tasks would clobber
each other if they shared `INP_OUT/`.

---

## 2. SLURM submission pattern

Production array job, amilan, MATLAB+libRadtran:

```bash
#SBATCH --account=ucb762_asc1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=82G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --array=1-1000%48
#SBATCH --chdir=/projects/anbu8374/slurm_logs/29_April_2026_synthetic
#SBATCH --output=create_meas_synthetic_NN_trainingData_BATCH_%A_%a.out
#SBATCH --error=create_meas_synthetic_NN_trainingData_BATCH_%A_%a.err
```

Things worth knowing:

- **`--output`/`--error` are resolved against the directory `sbatch` was
  invoked from, before any `cd` inside the script runs.** Without
  `--chdir`, the `.out`/`.err` files land in `~` and chew through
  home-quota inodes (we accumulated ~12,000 small log files across six
  batches before catching this). SLURM does **not** create `--chdir`'s
  target; `mkdir -p` it once before submitting.
- **Throttle concurrency**: `--array=1-1000%48` caps to 48 simultaneous
  tasks. Lower (e.g. `%40`) if Lustre contention reappears.
- **Chunk math** (we batch 57 clouds per task; 6 batches × 1000 tasks
  cover 300,001 clouds; last batch shrinks to `--array=1-248%48` with a
  `CLOUD_OFFSET` bump):
  ```
  start = CLOUD_OFFSET + (t-1)*CHUNK_SIZE + 1
  end   = min(CLOUD_OFFSET + t*CHUNK_SIZE, N_TOTAL)
  ```

Testing variant (interactive iteration):

```bash
#SBATCH --time=01:00:00
#SBATCH --partition=atesting
#SBATCH --qos=testing
#SBATCH --mem=50G
#SBATCH --cpus-per-task=10
#SBATCH --array=1
```
Note: `atesting` has been observed to grant more CPUs than asked for
(asked 10, got `SLURM_CPUS_PER_TASK=14`) — read the env var, don't
assume.

---

## 3. MATLAB-side parallelization

The driver loops over 57 clouds inside one MATLAB process; each cloud
calls `hysics_refl_from_synthetic_NN_inputs(...)` which uses `parfor`
internally to fan out `uvspec` invocations.

`start_parallel_pool.m` (curc branch) handles three things that
matter:

**1. Pool sized from `SLURM_CPUS_PER_TASK`, not `parcluster.NumWorkers`.**
`parcluster('local').NumWorkers` returns the **physical core count of the
node** (e.g. 64 on amilan), not the SLURM allocation. Sizing the pool
from that oversubscribes the cgroup and tanks throughput. Read the env
var instead:

```matlab
slurm_cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isnan(slurm_cpus) && slurm_cpus > 0
    num_workers = slurm_cpus;
end
parpool(p, num_workers);          % pass cluster object p, not just an int
```

**2. `JobStorageLocation` redirected to `$TMPDIR` (node-local).**
The default is `~/.matlab/local_cluster_jobs/<release>/`, which sits on
home (slow Lustre + tiny quota) and is shared across concurrent array
tasks, which corrupts pool state:

```matlab
tmpdir = getenv('TMPDIR');
if ~isempty(tmpdir)
    p.JobStorageLocation = tmpdir;
end
```
The bash script sets `$TMPDIR` to `/tmp/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}`
before MATLAB launches.

**3. Pool reuse + longer startup timeout.**
Tearing the pool down and recreating it costs minutes on amilan, so the
function returns early if a healthy pool of the right size already
exists. We also extend the start timeout from 20 → 40 min:

```matlab
pctconfig('poolstarttimeout', seconds(2400));
```

Observed pool startup once `$TMPDIR` is node-local: **30–60 s**.
Previously on Lustre, two tasks (481, 525) stalled multiple minutes in
"Job Queued" before workers attached.

---

## 4. MATLAB-specific issues we hit

**a. `mkdir: cannot create directory '/home/anbu8374/.MathWorks': Not a directory`** — appears in every `.err`. CURC explicitly disables this directory in their MATLAB wrapper; their own banner says to ignore. It's harmless.

**b. Multi-line `matlab -r "..."` strings break under CURC's MATLAB module wrapper.** The wrapper doesn't preserve embedded newlines cleanly and MATLAB ends up parsing a truncated `-r` argument. **Fix: collapse to one long line with `;` separators.** The pattern that works in `synthetic_4_May_2026.sh` looks like:
```bash
time matlab -nodesktop -nodisplay -r "addpath(genpath('...')); ...; start_parallel_pool('curc'); for cloud_id = ${start_id}:${end_id}, try, hysics_refl_from_synthetic_NN_inputs(...); catch ME, fprintf('...', ME.message); end; end; exit"
```

**c. Home directory bloat from `~/.matlab/local_cluster_jobs/`.** Grew
to 1.4 GB. Even with the script overriding `$TMPDIR`, any
interactive/ad-hoc parpool still falls back to home. Belt-and-suspenders
fix:

```bash
export MATLAB_PREFDIR=/projects/$USER/.matlab_prefs/R2024b
export MATLAB_LOG_DIR=/projects/$USER/.matlab/logs
mkdir -p "$MATLAB_PREFDIR" "$MATLAB_LOG_DIR"
```

```matlab
c = parcluster('local');
c.JobStorageLocation = '/projects/anbu8374/.matlab_prefs/local_cluster_jobs/R2024b';
saveProfile(c);
```
Plus a `find ... -mtime +7 -delete` prune in the trailer.

**d. `seff` is broken on `login-ci5`.** `Can't locate Slurmdb.pm in @INC` from `/curc/sw/slurmtools/latest/bin/seff`. Use a different login node (`login-ci3` works) or `sacct` directly.

---

## 5. Bottlenecks identified and how we beat them

The dominant bottleneck was **Lustre metadata-server contention**, not
CPU and not bandwidth. libRadtran's I/O pattern — open small file →
write a few KB → close → spawn `uvspec` → read it back → unlink — is
the worst-case workload for Lustre's MDS. At 48 array tasks × 40 parfor
workers × hundreds of clouds, every `open`/`stat`/`unlink` hits the
single shared MDS.

**Fixes, in the order they paid off:**

1. **Move `INP_OUT/`, `wc/`, `atmmod/` to node-local `/tmp`.** `/tmp` on
   amilan is a real ~63 GB SSD per node (confirmed via `df -h`), not
   tmpfs. The script tries `/tmp` first and falls back to `/scratch` only
   if the node-local mkdir fails:
   ```bash
   TMPDIR_LOCAL="/tmp/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
   if mkdir -p "$TMPDIR_LOCAL" 2>/dev/null && [ -w "$TMPDIR_LOCAL" ]; then
       export TMPDIR="$TMPDIR_LOCAL"
   else
       export TMPDIR="/scratch/alpine/${USER}/matlab_tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
       mkdir -p "$TMPDIR"
   fi
   ```
   A single `rm -rf "$TMPDIR"` in the trailer cleans up everything
   (INP/OUT, wc, atmmod, and the parpool JobStorageLocation files).

2. **Keep Mie tables on `/scratch`.** They're write-once, read-many
   across the whole campaign; putting them on `/tmp` would force
   regeneration on every task. Pre-create the parent once:
   ```bash
   mkdir -p "/scratch/alpine/${USER}/Mie_Calculations/"
   ```

3. **Stagger MATLAB starts** to soften the module-load thundering herd
   on `/projects`:
   ```bash
   sleep $((SLURM_ARRAY_TASK_ID % 10))
   ```

4. **Cap simultaneous array tasks.** `%48` was chosen empirically as a
   middle ground; the prior `%50` showed mild scratch-FS contention.

5. **Periodic prune** of stale node-local fallback dirs and per-task
   scratch dirs (7-day mtime cutoff) so a failed cleanup never
   accumulates.

**Throughput numbers (job 26685405, the first full production batch):**

| Configuration                                  | Per-cloud time | Notes |
|------------------------------------------------|----------------|-------|
| atesting, 14 workers, Lustre INP_OUT           | ~415 s         | hit wall before all 57 clouds finished |
| amilan, 40 workers, **node-local** `$TMPDIR`   | **~180 s**     | full 57 clouds in 142–180 min |

Speedup from 14 → 40 workers was 2.31× (linear would be 2.86×) —
slightly sublinear, consistent with per-cloud serial overhead.

Batch 26685405 final: **57,000 clouds ok / 0 failed across 1000 tasks**,
wallclock per task ranging 142–708 min (the slow tail is real but well
inside the 24 h budget).

---

## 6. Miscellaneous gotchas worth remembering

- **`/tmp` is per-node.** Two array tasks see two different `/tmp`s.
  Files there are guaranteed to survive until the task exits; lifetime
  after that is site-policy-dependent and we haven't pinned the exact
  number.
- **`curc-quota` is cached.** Freshly freed space won't show for several
  minutes; trust `du -sh` instead.
- **al40 (GPU) jobs occasionally land in "launch failed requeued
  held"** state — cause not yet diagnosed. Affects the NN sweep, not
  libRadtran itself.
- **conda env on al40**: moving the env from `~` to `/projects` broke
  `conda activate <name>`. Use the full path:
  `conda activate /projects/$USER/software/anaconda/envs/dropProfs_nn`.
- **Output filename convention** from the MATLAB driver encodes
  geometry, e.g.
  `simulated_spectra_HySICS_reflectance_636bands_0.3%_uncert_syntheticORACLES_cloud909_sza26_saz26_vza6_vaz56_sim-ran-on-04-May-2026.mat`.

---

## Key files

- `Hyperspectral_Cloud_Retrievals/Batch_Scripts/bash_scripts/neural_network/training_data/29_April_2026_synthetic/synthetic_4_May_2026.sh` — production array script
- `Hyperspectral_Cloud_Retrievals/Batch_Scripts/bash_scripts/.../synthetic_4_May_2026_TEST.sh` — atesting variant
- `startup/start_parallel_pool.m` — pool sizing + `JobStorageLocation` override
- `Hyperspectral_Cloud_Retrievals/HySICS/define_folderPaths_for_HySICS.m` — path map, `curc` branch
