# Daftar gambar yang disarankan untuk pembahasan pengaruh pH terhadap kristalisasi silika

## A. Gambar eksperimen dari skripsi yang perlu ditampilkan ulang/diringkas

| No. | Gambar | Peran dalam narasi |
|---:|---|---|
| 1 | Spektra FTIR variasi pH sebelum/atau setelah pengeringan | Menunjukkan evolusi Si-O-Si, Si-OH, dan air terikat. Narasi utama: kenaikan pH memperkuat kondensasi siloksan, tetapi basa berlebih dapat mempercepat aglomerasi. |
| 2 | XRD sebelum kalsinasi untuk pH 6,0; 6,5; 7,0; 7,5; 8,0 | Bukti bahwa seluruh prekursor masih dominan amorf, sehingga kristalisasi terjadi terutama selama kalsinasi/annealing. |
| 3 | Raman xerogel menurut pH | Menjelaskan cacat cincin kecil, keteraturan lokal, dan kontribusi gugus silanol pada prekursor amorf. |
| 4 | DTA/TG pH 6 | Menggambarkan prekursor sisi asam: air/gugus volatil lebih mudah lepas dan stabilitas termal lebih rendah. |
| 5 | DTA/TG pH 7 | Menjadi acuan kondisi optimum: dekomposisi lebih terkendali dan reorganisasi termal lebih stabil. |
| 6 | DTA/TG pH 8 | Menunjukkan sisi basa: stabilitas air terikat lebih rendah/tinggi tergantung tahap, serta indikasi reorganisasi fase pada suhu tinggi. |
| 7 | XRD setelah kalsinasi 900 °C | Gambar kunci: kemunculan cristobalite dan tridymite serta intensitas tertinggi pada pH optimum. |
| 8 | Plot turunan dari Tabel 4.1: pH vs ukuran kristalit cristobalite/tridymite | Sudah dibuat di `outputs/thesis_reference_plots/fig_thesis_crystallite_size_vs_pH.png`. Narasi: ukuran maksimum pada pH 7. |
| 9 | Plot turunan dari Tabel 4.1: pH vs rasio Ik/It | Sudah dibuat di `outputs/thesis_reference_plots/fig_thesis_Ikratio_vs_pH.png`. Narasi: pH 7 memberi dominasi cristobalite relatif paling kuat. |
| 10 | Plot turunan dari Tabel 4.1: pH vs FWHM | Sudah dibuat di `outputs/thesis_reference_plots/fig_thesis_fwhm_vs_pH.png`. Narasi: FWHM minimum pada pH 7 menandakan keteraturan kristal lebih baik. |

## B. Gambar simulasi yang perlu dibuat dari LAMMPS/ReaxFF/PLUMED

| No. | Gambar | File/analisis yang disiapkan | Makna pembahasan |
|---:|---|---|---|
| 11 | Skema workflow eksperimen-komputasi | Buat manual dari diagram alir README | Menghubungkan sol-gel, pH, kalsinasi, ReaxFF, XRD simulasi, dan MetaD. |
| 12 | Komposisi awal pH proxy: jumlah H2O/OH awal per pH | `structures/precursor_manifest_N1536.csv` | Menegaskan bahwa pH dimodelkan sebagai kondisi prekursor terhidrasi/terdefek, bukan asam-basa ionik eksplisit. |
| 13 | Temperatur, energi, tekanan, dan densitas selama annealing | Dari `log.lammps` setiap stage | Memastikan sistem stabil dan protokol termal tidak menghasilkan artefak. |
| 14 | Fraksi Si tetrahedral SiO4 terhadap pH | `scripts/05_structural_metrics.py` dan `08_plot_summary.py` | Mengukur pemulihan/kerusakan jaringan tetrahedral selama kalsinasi. |
| 15 | Spesiasi oksigen: bridging O, non-bridging O, silanol, water proxy | `fig_md_oxygen_speciation_vs_pH.png` | Menjelaskan mekanisme: asam menyisakan silanol, netral memaksimalkan siloksan, basa kuat dapat mempercepat pembentukan defek/agregasi. |
| 16 | Overlay Debye-XRD simulasi antar pH | `fig_md_xrd_overlay_vs_pH.png` | Membandingkan puncak/halo MD dengan XRD eksperimen. |
| 17 | Intensitas puncak simulasi sekitar 2θ ≈ 21-22° dan 35-36° vs pH | Turunan dari `04_compute_debye_xrd.py` | Menghubungkan order parameter simulasi dengan cristobalite-like ordering. |
| 18 | Snapshot atomistik prekursor dan struktur pasca-kalsinasi untuk pH 6, 7, 8 | OVITO/VMD dari `*.lammpstrj` | Visualisasi naratif: pH 7 lebih terkondensasi dan homogen. |
| 19 | Distribusi koordinasi Si-O dan O-H sebelum/sesudah kalsinasi | `05_structural_metrics.py` atau analisis tambahan | Bukti kuantitatif dehidroksilasi dan pembentukan Si-O-Si. |
| 20 | Local-Q6/CV trajectory dan HILLS MetaD | Stage 05, file `COLVAR` dan `HILLS_*` | Menunjukkan apakah enhanced sampling benar-benar mengeksplorasi keteraturan kristalin. |
| 21 | Free Energy Surface terhadap CV kristalisasi | PLUMED `sum_hills`/reweighting | Gambar puncak pembahasan: pH optimum diharapkan memiliki barrier kristalisasi lebih rendah. |
| 22 | Ring-size distribution jaringan Si-O | R.I.N.G.S./OVITO/skrip tambahan | Menguji perubahan topologi cincin menuju motif cristobalite/tridymite. |

## C. Urutan narasi pembahasan yang paling kuat

1. Mulai dari bukti eksperimen: semua prekursor amorf sebelum kalsinasi, tetapi pasca-kalsinasi menghasilkan cristobalite/tridymite.
2. Tunjukkan bahwa pH 7 adalah optimum eksperimen: ukuran kristalit maksimum, FWHM minimum, dan Ik/It tertinggi.
3. Gunakan ReaxFF untuk menjelaskan mekanismenya: keseimbangan silanol-siloksan, koordinasi SiO4, dan berkurangnya defek jaringan.
4. Gunakan Debye-XRD simulasi untuk menutup gap antara koordinat atomistik dan pola XRD eksperimen.
5. Jika custom XRD-CV PLUMED tersedia, gunakan FES untuk menyatakan bahwa prekursor pH 7 memiliki jalur kristalisasi yang lebih menguntungkan secara termodinamik/kinetik.
