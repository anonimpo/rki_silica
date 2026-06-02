# Catatan urutan definisi LAMMPS-ReaxFF-PLUMED

Urutan yang dipakai dalam semua input LAMMPS di direktori `lammps/` adalah:

1. `units real`, `atom_style charge`, dan `boundary p p p`.
2. `read_data` untuk membaca tipe atom dan koordinat.
3. `pair_style reaxff ...`.
4. `pair_coeff * * potentials/ffield.reax.SiOH Si O H`.
5. `neighbor` dan `neigh_modify`.
6. `fix qeq/reaxff ...` sebelum `run` atau `minimize`.
7. `timestep`, `thermo`, `dump`, lalu integrator (`fix npt`/`fix nvt`) atau `minimize`.
8. Jika memakai PLUMED: `fix plumed` didefinisikan setelah gaya ReaxFF dan QEq sudah aktif, serta sebelum `run` produksi yang ingin dibias.

Pemetaan tipe atom di seluruh data file:

| Tipe LAMMPS | Unsur | Keterangan |
|---:|---|---|
| 1 | Si | atom silikon jaringan silika |
| 2 | O  | oksigen jaringan + air/hidroksil |
| 3 | H  | hidrogen air/hidroksil |

Konsekuensi penting: `pair_coeff` **tidak boleh** ditulis `H O Si` untuk data ini, karena urutan setelah nama file harus mengikuti tipe atom LAMMPS, bukan urutan unsur di dalam file force field. Untuk data ini yang benar adalah:

```lammps
pair_coeff * * potentials/ffield.reax.SiOH Si O H
```

