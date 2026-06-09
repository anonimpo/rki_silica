import torch
import ase.io
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics

print("1. Membaca file struktur data...")
atoms = ase.io.read("liquid_silica_2400K.data", format="lammps-data")

# --- MODIFIKASI DISINI UNTUK GPU ---
# Mengecek apakah GPU NVIDIA (CUDA) tersedia di komputer Anda
if torch.cuda.is_available():
    device = "cuda"
    print("-> GPU terdeteksi! Simulasi akan berjalan di GPU.")
else:
    device = "cpu"
    print("-> GPU tidak terdeteksi. Simulasi berjalan di CPU (lambat).")

print("2. Memuat model pretrained CHGNet ke " + device)
# Memuat model langsung ke perangkat (GPU/CPU) yang dipilih
chgnet = CHGNet.load(use_device=device)
# -----------------------------------

print("3. Mengonfigurasi simulasi Molecular Dynamics (NVT @ 300K)...")
md = MolecularDynamics(
    atoms=atoms,
    model=chgnet,
    ensemble="nvt",
    temperature=300,
    timestep=0.25,
    trajectory="md_out.traj",
    logfile="md_out.log",
    loginterval=4,
)

print("4. Menjalankan simulasi MD... Silakan tunggu.")
md.run(40000)
print("   Simulasi MD Berhasil!")

print("5. Mengonversi hasil simulasi ke format XYZ...")
traj_configs = ase.io.read("md_out.traj", index=":")
ase.io.write("trajectory_ftir.xyz", traj_configs, format="extxyz")
