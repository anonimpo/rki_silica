# Asumsi model pH dalam paket ini

Force field yang diberikan (`ffield.reax.SiOH`) hanya memuat unsur H/O/Si. Karena itu, paket ini **belum** memasukkan spesies eksplisit seperti Na⁺, Cl⁻, HCl, NaOH, atau counter-ion lain. Representasi pH yang dipakai adalah **pH-conditioned precursor model**:

- pH 6,0 dan pH 8,0 diberi muatan hidrasi/defek awal lebih tinggi;
- pH 6,5 dan pH 7,5 berada pada tingkat antara;
- pH 7,0 diberi tingkat hidrasi/defek paling rendah, sesuai tren eksperimen bahwa pH 7 menghasilkan ukuran kristalit dan rasio Ik/It tertinggi setelah kalsinasi.

Makna ilmiah model ini adalah menguji bagaimana kondisi awal jaringan silika yang lebih/kurang terhidroksilasi memengaruhi densifikasi, dehidroksilasi, pembentukan jembatan Si-O-Si, dan kecenderungan keteraturan kristalin.

Untuk menguji **asam-basa eksplisit**, diperlukan salah satu rancangan berikut:

1. ReaxFF yang tervalidasi untuk Si/O/H/Na/Cl jika pH diatur melalui NaOH/HCl.
2. Model hidronium/hidroksida eksplisit dengan netralitas muatan total dan protokol jumlah ion yang konsisten dengan volume simulasi.
3. Validasi tambahan terhadap data eksperimen: XRD mentah, FTIR/Raman mentah, massa DTA/TG, serta komposisi kimia residu.

Dengan ukuran box MD nanometrik, istilah pH bulk tidak dapat diterapkan secara langsung tanpa kehati-hatian karena jumlah ion yang merepresentasikan pH tertentu sering kali menjadi fraksional. Oleh sebab itu, pH di sini diperlakukan sebagai variabel preparatif yang mengubah kimia prekursor, bukan sebagai larutan bulk ideal.
