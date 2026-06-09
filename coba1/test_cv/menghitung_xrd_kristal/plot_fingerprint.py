import numpy as np
import matplotlib.pyplot as plt

# 1. Tentukan nama file data fingerprint PLUMED Anda
filename = 'data_fingerprint_kristal.txt'

try:
    # 2. Baca data teks
    # np.loadtxt secara otomatis mengabaikan baris yang diawali dengan '#' atau '#!'
    data = np.loadtxt(filename)
    
    # Ekstraksi kolom sesuai format FIELDS PLUMED
    # Kolom 0 = q6_avg (Sumbu X)
    # Kolom 1 = hh2 / Histogram Densitas (Sumbu Y)
    q6_avg = data[:, 0]
    density = data[:, 1]

    # 3. Membuat plot grafik
    plt.figure(figsize=(9, 5.5))
    
    # Plot kurva fingerprint
    # Menggunakan fill_between untuk memberi warna transparan di bawah kurva agar lebih estetis
    plt.plot(q6_avg, density, color='#2c3e50', linewidth=2, label='Fase Kristal (Cristobalite)')
    plt.fill_between(q6_avg, density, color='#3498db', alpha=0.3)
    
    # 4. Mempercantik tampilan grafik
    plt.title('Fingerprint Struktural Silika - Distribusi Parameter $Q_6$', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Parameter Order Steinhardt Lokal ($Q_6$ Avg)', fontsize=12)
    plt.ylabel('Kerapatan Probabilitas (Density)', fontsize=12)
    
    # Rentang batas sumbu X ditentukan dari min_q6_avg (0.0) sampai max_q6_avg (1.0) sesuai header file Anda
    plt.xlim(0.0, 1.0)
    
    # Mengatur agar batas bawah sumbu Y pas di angka 0
    plt.ylim(bottom=0)

    # Menambahkan grid/garis bantu transparan
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(fontsize=11, loc='upper right')
    plt.tight_layout()

    # Menyimpan hasil plot ke dalam gambar PNG kualitas HD
    output_image = 'grafik_fingerprint_q6.png'
    plt.savefig(output_image, dpi=300)
    print(f"Sukses! Grafik berhasil disimpan sebagai '{output_image}'")

    # Menampilkan grafik di layar
    plt.show()

except FileNotFoundError:
    print(f"Error: File '{filename}' tidak ditemukan. Pastikan file berada di folder yang sama dengan skrip python ini.")
except Exception as e:
    print(f"Terjadi kesalahan saat membaca atau memplot data: {e}")
