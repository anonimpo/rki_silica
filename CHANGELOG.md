# Changelog

All notable changes to this project are documented below.

## [2026-06-07]

### Added
- **Pause Capability**: Added a graceful pipeline termination mechanism by checking for the presence of a `pause_flag` file at the start of [02_run_one_pipeline.sh](file:///home/rfa/repo/rki/rki_silica/scripts/02_run_one_pipeline.sh). Creating this file aborts pipeline execution with exit code 99 before starting the next pH case.
- **Aggregated N192 Output**: Generated XRD overlays and structural speciation plots for the $N=192$ system (`pH6p0`, `pH6p5`, `pH7p0`) under `/outputs/03_calcination_1173K_N192`.

### Changed
- **Parallel Optimization**: Unquoted the `$LMP` variable expansion in [02_run_one_pipeline.sh](file:///home/rfa/repo/rki/rki_silica/scripts/02_run_one_pipeline.sh) to support multi-word execution commands (e.g., `LMP="mpirun -np 16 lmp"`). This enables MPI-based spatial domain decomposition which yields a >2x speedup on small systems and near-linear scaling on larger ones.
- **Python Dependencies**: Configured missing packages (`pandas` and `matplotlib`) in the conda environment to resolve post-processing script import failures.

### Ongoing Work
- **$N=1536$ Production Run**: Successfully completed Stage 03 calcination for `pH6p0` on 16 cores. Resumed Stage 04 high-temperature production manually. The run is configured to pause after `pH6p0` is fully finished.

### How to Resume/Continue Calculations

1. **Resume Stage 04 for `pH6p0`**:
   If Stage 04 was paused/cancelled, resume it using:
   ```bash
   mpirun -np 16 lmp -var phlabel pH6p0 -var baseatoms 1536 -var Tprod 2300 -var nheat 80000 -var nprod 400000 -in lammps/in.04_highT_unbiased_reaxff
   ```

2. **Continue Remaining Cases (`pH6p5`, `pH7p0`, `pH7p5`, `pH8p0`)**:
   Once `pH6p0` is fully complete, delete the `pause_flag` file and run the remaining cases in the queue:
   ```bash
   # Remove the pause flag
   rm pause_flag

   # Run remaining cases sequentially
   for PH in pH6p5 pH7p0 pH7p5 pH8p0; do
     LMP="mpirun -np 16 lmp" ./scripts/02_run_one_pipeline.sh "$PH" 1536
   done
   ```
