# silica_pH7_N192_20core

Direktori kerja ini dibuat untuk percobaan pH 7, base silica 192 atom dengan LAMMPS + ReaxFF dan opsi PLUMED. Pada file struktur, base SiO2 berisi 64 Si + 128 O = 192 atom. Karena model pH 7 ditambah 8 molekul H2O, total atom aktual dalam file LAMMPS menjadi 216 atom.

Struktur folder:

- `kode_lammps/` berisi input LAMMPS, file PLUMED, struktur awal, dan force field.
- `post_processing/` berisi skrip XRD Debye dan visualisasi snapshot mirip Fig. 7 Niu et al.
- `output/` berisi hasil plot.
- `runs/` berisi hasil run. Saya sertakan baseline pH7_N192 dari pekerjaan sebelumnya supaya post-processing bisa langsung diuji.

## Cara menjalankan dengan 20 core

Dari root folder proyek:

```bash
chmod +x kode_lammps/run_pH7_N192_20core.sh
NP=20 LMP=lmp ./kode_lammps/run_pH7_N192_20core.sh
```

Jika executable LAMMPS bernama lain:

```bash
NP=20 LMP=lmp_mpi ./kode_lammps/run_pH7_N192_20core.sh
```

## Catatan penting PLUMED

Error lama disebabkan urutan fix berikut:

```lammps
fix prod all npt ...
fix plm all plumed ...
```

Pada `kode_lammps/in.05_wtmetad_plumed_FIXED`, urutannya sudah diperbaiki menjadi:

```lammps
fix plm all plumed plumedfile ${plumedfile} outfile ${rundir}/plumed.log
fix prod all npt ...
```

Dengan urutan ini, `fix plumed` didefinisikan sebelum `fix npt` pada tahap produksi.

## Cara membuat output gambar

Dari root folder proyek:

```bash
pip install -r post_processing/requirements.txt
bash post_processing/run_postprocess.sh
```

Output utama:

- `output/xrd_before_after_pH7_N192.png`
- `output/figure7_like_pH7_N192.png`

## Batasan ilmiah

Model 192 atom adalah percobaan kecil untuk debugging dan validasi alur kerja. Untuk klaim kristalisasi seperti Niu et al., sistem kecil dan waktu simulasi pendek biasanya belum cukup. Niu et al. menggunakan sistem jauh lebih besar dan enhanced sampling. Karena itu, output dari folder ini sebaiknya diperlakukan sebagai workflow awal, bukan hasil final paper.
