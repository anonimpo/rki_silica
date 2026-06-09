# Niu et al. WTMETAD silica crystallization reproduction workspace

Target:
1. make/check XRD before crystallization: liquid/amorphous-like broad pattern;
2. run WTMETAD at 2300 K with an XRD {111}-like CV;
3. make XRD after WTMETAD and compare whether sharp beta-cristobalite-like peaks appear;
4. make a Fig. 7-like snapshot panel with only crystal-like Si atoms shown.

Default system size follows Niu Fig. 7: **1536 atoms = 512 Si + 1024 O**.
Default run uses **20 MPI cores**.

## Important honesty note

The exact Niu paper used an XRD peak-intensity CV based on the Debye scattering equation with atomic scattering form factors. Standard PLUMED 2.10 documentation does not show a ready-made `XRD_INTENSITY` action. Therefore this workspace provides:

- `custom_cv/plumed_xrd111_si_matrix_WTMETAD_v210.dat`: a documented-PLUMED-2.10, matrix-based **Si-only Debye-style XRD {111} CV**. This is the closest reproducible PLUMED-native test included here.
- `custom_cv/plumed_envsimilarity_fallback_WTMETAD_v210.dat`: fallback using PLUMED-native `ENVIRONMENTSIMILARITY`. This is **not** Niu's exact CV.
- `postprocessing/compute_xrd_debye.py`: postprocessing Debye XRD for before/after comparison.

If the research target is a strict one-to-one reproduction, the next step is to implement a true custom PLUMED action with all-atom XRD form factors and analytic derivatives.

## Quick order

```bash
unzip niu_wtmetad_repro_20core.zip
cd niu_wtmetad_repro_20core
chmod +x runs/*.sh postprocessing/*.sh postprocessing/*.py

# check installation
LMP=lmp NCORE=20 ./runs/00_check_environment.sh

# parse test for LAMMPS + PLUMED + XRD-like CV
LMP=lmp NCORE=20 ./runs/02_smoke_test_plumed_xrd_20core.sh

# make liquid/amorphous-like starting structure from beta-cristobalite
LMP=lmp NCORE=20 ./runs/03_make_liquid_from_beta_20core.sh

# optional unbiased control
LMP=lmp NCORE=20 ./runs/04_run_unbiased_control_20core.sh

# main WTMETAD check
LMP=lmp NCORE=20 ./runs/05_run_wtmetad_xrd111_20core.sh

# XRD before/after + Fig. 7-like image
./runs/07_postprocess.sh
```

For quick functional testing only:

```bash
LMP=lmp NCORE=20 ./runs/run_all_quick_check_20core.sh
```

The quick check is intentionally short and may not crystallize. It only checks whether the workflow runs.

## Expected output files

- `outputs/COLVAR_XRD111_WTMETAD`: time series of the XRD-like CV and metadynamics bias.
- `outputs/HILLS_XRD111`: Gaussian hills from WTMETAD.
- `outputs/xrd_before_after_wtmetad_check.png`: before/after XRD comparison.
- `outputs/fig7_like_wtmetad_xrd111.png`: Fig. 7-like six-panel visualization if the WTMETAD dump exists.
- `structures/wtmetad_xrd111_after_2300K_N1536.data`: final structure after WTMETAD.

## Default physical choices

Niu used a 2 fs timestep with the Takada silica potential. This workspace uses ReaxFF because your setup has ReaxFF, and the default timestep is **0.25 fs** for stability. Therefore Gaussian deposition every 1 ps is set as `PACE=4000` in PLUMED. If you intentionally use 2 fs, change `PACE=500` in `custom_cv/plumed_xrd111_si_matrix_WTMETAD_v210.dat`.

## Known failure modes

1. `Unknown fix style plumed`: LAMMPS was not built with the PLUMED package.
2. `DISTANCE_MATRIX` not recognized: PLUMED was not configured with the `adjmat` module.
3. `PAIRENTROPIES` or `LOCAL_AVERAGE` syntax error: use the v2.10-adapted file, not the original SI syntax.
4. ReaxFF instability at high T: reduce timestep, reduce melt temperature, or increase damping.
5. No sharp peak after short run: increase `RUN_STEPS`; crystallization is still a rare event and the quick run is not enough for scientific conclusion.
