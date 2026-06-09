# Documentation audit and non-hallucination notes

Verified from Niu et al. paper/SI:

- The paper uses XRD peak intensity as collective variable, especially the beta-cristobalite `{111}` peak.
- The Debye scattering expression is used for powder-like XRD intensity.
- The paper's computational setup reports: LAMMPS + PLUMED 2 development version, 2 fs timestep, WTMETAD bias factor 100, Gaussian deposition every 1 ps, width 5 CV units, height 40 kJ/mol.
- Fig. 7 uses 1536 atoms total = 512 Si atoms, 2300 K, and snapshots every 40 ps.
- The SI local entropy fingerprint uses only Si atoms and parameters rm=0.75 nm, ra=0.45 nm, sigma=0.05 nm.

Verified from PLUMED 2.10 documentation:

- `METAD` is the documented action for metadynamics and WTMETAD; it writes Gaussian history to `HILLS` and uses keywords such as `ARG`, `SIGMA`, `HEIGHT`, `PACE`, `BIASFACTOR`, `TEMP`, `GRID_MIN`, `GRID_MAX`, `GRID_BIN`.
- `LOAD` can load a shared object or a C++ source file defining new actions. This is the correct path for a future exact custom XRD CV action.
- `DISTANCE_MATRIX` is part of the `adjmat` module and can compute distance/component matrices.
- `CUSTOM` can combine variables through algebraic expressions and automatically differentiates them.
- `SUM` computes sums of arguments.
- `ENVIRONMENTSIMILARITY` is part of the `envsim` module and can compare local environments to a reference crystal structure, including DIAMOND.
- LAMMPS `fix plumed` syntax is `fix ID group-ID plumed plumedfile <file> outfile <file>`, and it should appear after relevant input parameters such as timestep.

What is not claimed:

- This package does **not** claim that standard PLUMED 2.10 has a built-in exact Niu XRD CV action.
- The included XRD-like WTMETAD CV is Si-only and omits atomic X-ray form factors. It is a practical test of the core idea, not a perfect reproduction.
- The fallback `ENVIRONMENTSIMILARITY` route is not the paper's method.
