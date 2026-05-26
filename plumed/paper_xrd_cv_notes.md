# Catatan CV XRD sesuai paper Niu et al. (2018)

Paper menggunakan collective variable berbasis intensitas puncak XRD:

- s1 = I(Q{111}) untuk bias/metadynamics.
- s2 = I(Q{022}) untuk analisis/reweighting.
- I(Q) dihitung sebagai Debye scattering function yang dinormalisasi jumlah atom dan diberi Lorch/window cutoff.
- Parameter WTMetaD: bias factor 100, Gaussian tiap 1 ps, sigma/width 5 unit CV, height 40 kJ/mol.

Penting: paper menyebut LAMMPS yang dipatch dengan development version PLUMED 2. Pada PLUMED stock v2.10, tidak ada action resmi bernama XRD/STRUCTUREFACTOR dalam indeks action; yang paling dekat adalah `SAXS`, tetapi itu bukan substitusi persis untuk CV XRD paper. Karena itu, paket ini menyediakan dua jalur:

1. **Jalur produksi paper-level**: compile satu executable LAMMPS dengan REAXFF + PLUMED + CV XRD/structure-factor khusus, lalu ubah `plumed_saxs_metad_approx.dat` ke action CV XRD Anda.
2. **Jalur praktis dua executable**: jalankan ReaxFF MD dengan `in.04_reaxff_production_no_plumed`, lalu hitung I(Q111) dan I(Q022) dengan `scripts/xrd_debye_cv.py` dari dump trajectory. Jalur ini tidak memberi bias balik ke dinamika, tetapi membantu validasi dan analisis.

Untuk q dari sel kubik beta-cristobalite a=7.15 A:

- Q111 = 1.52207 A^-1
- Q022 = 2.48553 A^-1
