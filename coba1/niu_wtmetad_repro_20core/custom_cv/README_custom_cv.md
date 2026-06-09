# Custom CV notes

## Main WTMETAD CV

`plumed_xrd111_si_matrix_WTMETAD_v210.dat` builds a Si-only Debye-style peak intensity using documented PLUMED 2.10 actions:

1. `GROUP` selects Si atoms 1-512.
2. `DISTANCE_MATRIX` calculates Si-Si distances.
3. `CUSTOM` evaluates the Debye term `sin(Qr)/(Qr)` multiplied by a Lorch-style window.
4. `SUM` sums pair contributions.
5. `METAD` biases the resulting scalar.

This tests the Niu idea directly: increase beta-cristobalite `{111}`-like order by biasing an XRD-like intensity.

## Limitation

The paper CV includes all atoms and atomic X-ray scattering form factors. This PLUMED-native file uses Si-Si only. It is intentionally transparent rather than pretending to be exact.

## Exact future route

For exact reproduction, implement `XRD_INTENSITY` as a PLUMED C++ action and load it with:

```plumed
LOAD FILE=XRDIntensity.cpp
s1: XRD_INTENSITY SPECIES=1-1536 Q=1.522 RCUT=12.0 FORMFACTORS=Si,O
```

PLUMED 2.10 documentation supports loading custom actions via `LOAD`, including direct `.cpp` loading or compiled shared objects.
