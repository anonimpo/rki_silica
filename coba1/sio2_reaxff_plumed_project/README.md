# SiO2 crystallization workflow: LAMMPS ReaxFF + PLUMED

Paket ini menyiapkan struktur awal, input LAMMPS, input PLUMED, dan skrip analisis untuk memodelkan kristalisasi SiO2 menuju beta-cristobalite dengan inspirasi dari Niu et al. (2018). File force field `ffield.reax.SiOH` yang Anda unggah sudah disalin ke `potentials/`.

## Struktur direktori

```text
sio2_reaxff_plumed_project/
├── README.md
├── potentials/
│   └── ffield.reax.SiOH
├── structures/
│   ├── beta_cristobalite_192.data
│   ├── beta_cristobalite_192.pdb
│   └── beta_cristobalite_192.xyz
├── lammps/
│   ├── in.00_minimize_reaxff
│   ├── in.01_melt_reaxff
│   ├── in.02_equilibrate_liquid_reaxff
│   ├── in.03_reaxff_plumed_wtmetad
│   └── in.04_reaxff_production_no_plumed
├── plumed/
│   ├── paper_xrd_cv_notes.md
│   ├── plumed_saxs_metad_approx.dat
│   ├── plumed_saxs_print_only.dat
│   └── q_values_beta_cristobalite.txt
├── scripts/
│   ├── check_lammps_features.sh
│   ├── make_beta_cristobalite.py
│   ├── make_beta_cristobalite_ase.py
│   ├── postprocess_existing_dump.sh
│   ├── run_reaxff_pipeline.sh
│   ├── run_reaxff_postprocess_workflow.sh
│   ├── run_wtmetad.sh
│   └── xrd_debye_cv.py
├── outputs/
│   └── xrd_initial.dat
└── logs/
```

## Struktur awal

`structures/beta_cristobalite_192.data` berisi 192 atom: 64 Si + 128 O. Ini dibuat sebagai supercell 2x2x2 dari unit konvensional beta-cristobalite ideal, dengan parameter kisi `a = 7.15 Å`. Skala ini mengikuti ukuran kecil 192 atom yang dilaporkan pada paper untuk FES awal/ilustrasi.

Untuk membuat ulang struktur tanpa dependency eksternal:

```bash
cd sio2_reaxff_plumed_project
python3 scripts/make_beta_cristobalite.py --a 7.15 --reps 2 2 2 --out-prefix structures/beta_cristobalite_192
```

Alternatif jika ASE terinstal:

```bash
pip install ase
python3 scripts/make_beta_cristobalite_ase.py --a 7.15 --reps 2 2 2 --out-prefix structures/beta_cristobalite_192_ase
```

Untuk menambahkan sedikit gangguan posisi sebelum minimisasi:

```bash
python3 scripts/make_beta_cristobalite.py --a 7.15 --reps 2 2 2 --jitter 0.02 --out-prefix structures/beta_cristobalite_192_jitter
```

## Catatan penting tentang paper dan ReaxFF

Paper Niu et al. menggunakan LAMMPS + development PLUMED 2, timestep 2 fs, thermostat stochastic velocity rescaling dengan waktu relaksasi 0.1 ps, barostat Parrinello-Rahman pada 1 atm dengan waktu relaksasi 10 ps, dan WTMetaD dengan bias factor 100, Gaussian setiap 1 ps, sigma 5 unit CV, dan tinggi 40 kJ/mol.

Paket ini memakai ReaxFF karena Anda mengunggah `ffield.reax.SiOH`. Karena ReaxFF reactive dan QEq biasanya lebih sensitif, timestep default diset konservatif `0.25 fs`, bukan 2 fs. Setelah energi stabil, Anda boleh uji bertahap 0.5 fs atau lebih, tetapi jangan langsung memakai 2 fs untuk ReaxFF tanpa validasi energi/temperatur.

## Jalur A - satu executable LAMMPS dengan ReaxFF + PLUMED

Ini jalur yang diperlukan untuk metadynamics ter-bias. Satu executable LAMMPS harus mengenali `pair_style reaxff`, `fix qeq/reaxff`, dan `fix plumed`.

Cek executable:

```bash
cd sio2_reaxff_plumed_project
LAMMPS_BIN=/path/to/lmp ./scripts/check_lammps_features.sh
```

Jalankan preparasi:

```bash
LAMMPS_BIN=/path/to/lmp ./scripts/run_reaxff_pipeline.sh
```

Lalu jalankan WTMetaD approximate:

```bash
LAMMPS_BIN=/path/to/lmp ./scripts/run_wtmetad.sh
```

Input PLUMED yang disediakan (`plumed/plumed_saxs_metad_approx.dat`) adalah pendekatan berbasis stock `SAXS` PLUMED, bukan exact paper CV. Untuk reproduksi paper-level, Anda perlu PLUMED dengan CV XRD/structure factor development atau implementasi CV sendiri, lalu ganti isi file PLUMED tersebut.

## Jalur B - dua executable terpisah: ReaxFF saja + analisis CV offline

Jika Anda hanya punya executable ReaxFF dan executable/instalasi PLUMED terpisah, bias metadynamics tidak dapat diterapkan balik ke dinamika ReaxFF. Jalur praktisnya:

1. Jalankan ReaxFF untuk membuat trajectory.
2. Hitung CV XRD `I(Q111)` dan `I(Q022)` dari dump trajectory memakai Python.

Contoh:

```bash
cd sio2_reaxff_plumed_project
LAMMPS_BIN=/path/to/lmp_reaxff ./scripts/run_reaxff_pipeline.sh
LAMMPS_BIN=/path/to/lmp_reaxff ./scripts/run_reaxff_postprocess_workflow.sh
```

Atau untuk dump yang sudah ada:

```bash
./scripts/postprocess_existing_dump.sh outputs/production_no_plumed.lammpstrj outputs/xrd_cv_from_dump.dat
```

Format output:

```text
# step I_Q1.52207 I_Q2.48553
0 ... ...
1000 ... ...
```

## Urutan input LAMMPS

1. `lammps/in.00_minimize_reaxff`  
   Minimisasi beta-cristobalite ideal.

2. `lammps/in.01_melt_reaxff`  
   Melt pada 4000 K dengan NPT. Default pendek 10 ps; produksi perlu diperpanjang.

3. `lammps/in.02_equilibrate_liquid_reaxff`  
   Equilibrate liquid silica pada 2400 K dan 1 atm.

4. `lammps/in.03_reaxff_plumed_wtmetad`  
   ReaxFF + PLUMED WTMetaD. Perlu satu executable gabungan.

5. `lammps/in.04_reaxff_production_no_plumed`  
   Produksi ReaxFF tanpa PLUMED, untuk workflow dua executable dan post-processing.

## Parameter Q untuk CV

Dengan `a = 7.15 Å`:

```text
Q111 = 1.52207 Å^-1
Q022 = 2.48553 Å^-1
```

Nilai ini ada di `plumed/q_values_beta_cristobalite.txt`. Script Python `xrd_debye_cv.py` memakai nilai ini sebagai default.

## Cara build LAMMPS yang benar secara konsep

Anda butuh satu LAMMPS yang dibangun dengan paket REAXFF dan PLUMED. Contoh CMake konseptual:

```bash
git clone https://github.com/lammps/lammps.git
cd lammps
mkdir build && cd build
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=on \
  -D PKG_REAXFF=on \
  -D PKG_PLUMED=on
make -j 4
```

Jika PLUMED terinstal non-standar, pastikan environment PLUMED (`plumed config`) bisa ditemukan oleh CMake atau patch sesuai instruksi versi PLUMED/LAMMPS Anda.

## Troubleshooting cepat

- `ERROR: Unrecognized pair style 'reaxff'`: LAMMPS tidak punya paket REAXFF, atau binary yang dipakai salah.
- `ERROR: Unrecognized fix style 'plumed'`: LAMMPS tidak dipatch/dibuild dengan PLUMED.
- `Pair reaxff requires fix qeq/reaxff`: pastikan baris `fix qeq all qeq/reaxff ... reaxff` aktif.
- Simulasi meledak pada suhu tinggi: kecilkan timestep ke 0.1 fs, mulai dari minimisasi, atau panaskan bertahap.
- PLUMED gagal pada `ARG=xrd.q-0`: cek `logs/plumed_wtmetad.log` atau jalankan print-only untuk melihat nama komponen output `SAXS` di versi PLUMED Anda.

## Validasi yang sudah dijalankan di sandbox

Di sandbox ini tidak tersedia executable LAMMPS/PLUMED maupun modul ASE. Yang berhasil dijalankan:

```bash
python3 scripts/make_beta_cristobalite.py
python3 scripts/xrd_debye_cv.py structures/beta_cristobalite_192.xyz
```

Hasil awal tersimpan di `outputs/xrd_initial.dat`.
