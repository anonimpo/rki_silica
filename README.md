# pH-Conditioned Silica Crystallization: ReaxFF + LAMMPS + PLUMED Workspace

Workspace ini disiapkan untuk mengembangkan riset lanjutan dari hasil skripsi tentang variasi pH sintesis silika sekam padi, lalu diuji dengan dinamika molekular reaktif. Paket ini memuat struktur awal, input LAMMPS berurutan, template PLUMED, skrip analisis XRD-Debye, skrip metrik struktur, dan daftar gambar untuk pembahasan.

## 1. Struktur direktori

```text
ph_silica_reaxff_md_project/
├── potentials/                     # ffield.reax.SiOH dari file yang Anda unggah
├── structures/                      # prekursor pH 6,0–8,0 untuk N=192 dan N=1536
├── lammps/                          # input LAMMPS Stage 01–05
├── plumed/
│   ├── stock_localq6/               # PLUMED Local-Q6 fallback, stock PLUMED + crystallization module
│   └── custom_xrd_template/         # template XRD/Debye CV ala Niu; perlu custom action
├── scripts/                         # generator, runner, post-processing
├── outputs/thesis_reference_plots/   # plot dari Tabel 4.1 skripsi
├── docs/                            # catatan model, urutan definisi, dan daftar gambar
├── runs/                            # output LAMMPS akan masuk ke sini
└── analysis/                        # output post-processing akan masuk ke sini
```

## 2. Batasan penting model pH

Force field `ffield.reax.SiOH` hanya memuat H/O/Si. Karena itu, workspace ini belum memodelkan HCl, NaOH, Na⁺, atau Cl⁻ secara eksplisit. pH direpresentasikan sebagai **kondisi prekursor terhidrasi/terdefek**: pH 7,0 diberi defek/hidrasi awal terendah, sedangkan pH 6,0 dan 8,0 lebih tinggi. Rancangan ini selaras dengan tren skripsi bahwa pH 7 menghasilkan kristalit terbesar dan FWHM terendah setelah kalsinasi.

Untuk asam-basa eksplisit, diperlukan force field tambahan yang tervalidasi untuk ion terkait, atau protokol hidronium/hidroksida yang menjaga netralitas total sistem.

## 3. Cek lingkungan

```bash
cd ph_silica_reaxff_md_project
LMP=lmp ./scripts/00_check_environment.sh
```

Minimal diperlukan:

- LAMMPS dengan paket/fungsi `REAXFF` dan `fix qeq/reaxff`.
- Untuk Stage 05: LAMMPS terhubung dengan PLUMED dan PLUMED memiliki module yang sesuai.
- Python 3 dengan `numpy`, `pandas`, dan `matplotlib`.

## 4. Membuat ulang struktur awal

Struktur sudah disertakan. Untuk membuat ulang:

```bash
python3 scripts/01_make_ph_precursors.py --base-atoms 1536
```

Untuk uji cepat:

```bash
python3 scripts/01_make_ph_precursors.py --base-atoms 192 --scale-box 1.20
```

## 5. Menjalankan pipeline baseline

Uji cepat tanpa menjalankan LAMMPS:

```bash
DRYRUN=1 LMP=lmp ./scripts/03_run_all_ph.sh 192
```

Run satu pH:

```bash
LMP=lmp ./scripts/02_run_one_pipeline.sh pH7p0 1536
```

Run semua pH secara berurutan:

```bash
LMP=lmp ./scripts/03_run_all_ph.sh 1536
```

Untuk produksi ilmiah, naikkan jumlah step:

```bash
LMP=lmp \
N_EQ=200000 N_HEAT=400000 N_HOLD=800000 N_PROD=2000000 \
./scripts/03_run_all_ph.sh 1536
```

Dengan timestep 0,25 fs, 2.000.000 step = 500 ps.

## 6. Urutan stage LAMMPS

| Stage | File input | Tujuan |
|---:|---|---|
| 01 | `lammps/in.01_minimize_precursor` | minimisasi energi awal |
| 02 | `lammps/in.02_equilibrate_300K` | ekuilibrasi 300 K |
| 03 | `lammps/in.03_calcination_1173K` | analog kalsinasi 900 °C/1173 K |
| 04 | `lammps/in.04_highT_unbiased_reaxff` | trajektori suhu tinggi tanpa bias untuk screening kristalisasi dan rentang CV |
| 05 | `lammps/in.05_wtmetad_plumed` | WTMetaD opsional dengan PLUMED |

Semua input memakai urutan ReaxFF yang sama: `read_data` → `pair_style reaxff` → `pair_coeff * * ... Si O H` → `fix qeq/reaxff` → integrator/PLUMED → `run`.

## 7. PLUMED: dua opsi

### Opsi A — Stock Local-Q6 fallback

File tersedia di:

```text
plumed/stock_localq6/plumed_localq6_pH7p0_N1536.dat
```

Ini dapat dipakai untuk eksplorasi awal jika PLUMED dibangun dengan `crystallization` module. CV ini bukan CV XRD Niu, tetapi berguna untuk memantau peningkatan order pada sublattice Si.

Menjalankan Stage 05:

```bash
LMP=lmp RUN_META=1 ./scripts/02_run_one_pipeline.sh pH7p0 1536
```

### Opsi B — Niu-style XRD/Debye WTMetaD

Template ada di:

```text
plumed/custom_xrd_template/plumed_xrd_wtmetad_NEEDS_CUSTOM_DEBYE_ACTION.dat
```

Ini memerlukan custom PLUMED action `DEBYE_STRUCTURE_FACTOR` atau action ekuivalen. Jangan menjalankan file ini sebelum custom action benar-benar terkompilasi dan keyword-nya cocok dengan versi PLUMED yang dipakai.

## 8. Post-processing

Setelah Stage 03 selesai:

```bash
./scripts/07_postprocess_all.sh 1536 03_calcination_1173K calcined_1173K.data
```

Output utama:

```text
outputs/03_calcination_1173K_N1536/fig_md_xrd_overlay_vs_pH.png
outputs/03_calcination_1173K_N1536/fig_md_tetrahedral_si_fraction_vs_pH.png
outputs/03_calcination_1173K_N1536/fig_md_oxygen_speciation_vs_pH.png
analysis/03_calcination_1173K_N1536/metrics_all.csv
```

Plot referensi skripsi sudah dibuat di:

```text
outputs/thesis_reference_plots/
```

## 9. Cara membaca hasil untuk narasi pH → kristalisasi

Narasi inti yang disarankan:

1. Sebelum kalsinasi, eksperimen menunjukkan prekursor tetap amorf.
2. Setelah kalsinasi 900 °C, pH 7 menghasilkan kristalit cristobalite paling besar, FWHM terkecil, dan rasio Ik/It tertinggi.
3. ReaxFF digunakan untuk menguji apakah pH proxy tersebut menghasilkan jaringan dengan fraksi SiO4 lebih tinggi, silanol lebih rendah, dan jembatan Si-O-Si lebih dominan.
4. Debye-XRD simulasi digunakan untuk menghubungkan koordinat atomistik dengan puncak XRD eksperimen.
5. Jika custom XRD-CV WTMetaD tersedia, FES dapat dipakai untuk menguji apakah pH 7 memiliki barrier kristalisasi lebih rendah atau basin kristalin yang lebih stabil.

