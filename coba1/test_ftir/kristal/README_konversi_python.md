# Konversi notebook FTIR CHGNet MD ke Python

Isi file:

1. `run_md_chgnet.py`  
   Script utama untuk menjalankan MD ASE + CHGNet dan menyimpan `dipoles.csv`.

2. `ftir_postprocess.py`  
   Script untuk mengubah `dipoles.csv` menjadi spektrum FTIR relatif.

3. `ftir_chgnet_md_workflow.py`  
   Wrapper workflow agar konfigurasi dari notebook bisa dipanggil lewat command line.

## Instalasi dependency

```bash
pip install ase chgnet torch scipy matplotlib
```

Untuk GPU, gunakan versi `torch` yang sesuai dengan CUDA environment Anda.

## Contoh penggunaan

Cek input struktur:

```bash
python ftir_chgnet_md_workflow.py --inspect
```

Jalankan MD:

```bash
python ftir_chgnet_md_workflow.py --run-md --input /content/calcined_1173K.data --device cuda
```

Jalankan FTIR setelah `dipoles.csv` tersedia:

```bash
python ftir_chgnet_md_workflow.py --run-ftir
```

Uji post-processing FTIR tanpa MD:

```bash
python ftir_chgnet_md_workflow.py --self-test
```
