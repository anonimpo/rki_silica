import numpy as np
import matplotlib.pyplot as plt

# 1. Tentukan nama file data Anda
filename = 'spektrum_xrd.txt'

try:
    # 2. Baca data teks
    # skiprows=4 berfungsi untuk melewati 4 baris pertama (komentar & timestep)
    data = np.loadtxt(filename, skiprows=4)
    
    # Berdasarkan format LAMMPS:
    # Kolom 0 = Nomor Baris (Row)
    # Kolom 1 = Sudut 2-Theta (c_myXRD[1])
    # Kolom 2 = Intensitas (c_myXRD[2])
    
    two_theta = data[:, 1]
    intensity = data[:, 2]

    # Mengurutkan data berdasarkan nilai 2-Theta dari kecil ke besar 
    # (Penting agar garis grafik tidak zig-zag berantakan)
    sort_indices = np.argsort(two_theta)
    two_theta_sorted = two_theta[sort_indices]
    intensity_sorted = intensity[sort_indices]

    # 3. Membuat plot grafik
    plt.figure(figsize=(10, 6))
    
    # Anda bisa mengganti 'plot' menjadi 'scatter' jika ingin melihat titik-titik aslinya
    plt.plot(two_theta_sorted, intensity_sorted, color='blue', linewidth=1.5, label='Cristobalite')
    
    # 4. Mempercantik tampilan grafik
    plt.title('Simulasi X-Ray Diffraction (XRD) - Fase Kristal', fontsize=16, fontweight='bold')
    plt.xlabel(r'Sudut $2\theta$ (Derajat)', fontsize=14)
    plt.ylabel('Intensitas (a.u.)', fontsize=14)
    
    # Membatasi sumbu X sesuai input Anda (10 sampai 90 derajat)
    plt.xlim(10, 90)
    
    # Mengatur agar nilai Y (intensitas) terendah menempel di angka 0
    plt.ylim(bottom=0)

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    plt.tight_layout()

    # Menyimpan hasil plot ke dalam gambar PNG (opsional)
    plt.savefig('grafik_xrd.png', dpi=300)
    print("Grafik berhasil disimpan sebagai 'grafik_xrd.png'")

    # Menampilkan grafik di layar
    plt.show()

except FileNotFoundError:
    print(f"Error: File '{filename}' tidak ditemukan. Pastikan nama file dan foldernya benar.")
except Exception as e:
    print(f"Terjadi kesalahan: {e}")
