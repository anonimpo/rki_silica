import torch
import ase.io
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics, CHGNetCalculator  # <--- Tambahkan CHGNetCalculator
import os

print("1. Membaca file struktur data...")
data_file = "amorphous_pH7_N1536.data"
if not os.path.exists(data_file):
    potential_paths = [
        os.path.join("scripts", "custom", "amorphous_pH7_N1536.data"),
        os.path.join("rki_silica", "scripts", "custom", "amorphous_pH7_N1536.data"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "amorphous_pH7_N1536.data")
    ]
    for path in potential_paths:
        if os.path.exists(path):
            data_file = path
            break

print(f"-> Membaca dari: {data_file}")
atoms = ase.io.read(data_file, format="lammps-data")

# --- MODIFIKASI DISINI UNTUK GPU ---
if torch.cuda.is_available():
    device = "cuda"
    print("-> GPU terdeteksi! Simulasi akan berjalan di GPU.")
    # Mengurangi fragmentasi memori CUDA
    import os
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
else:
    device = "cpu"
    print("-> GPU tidak terdeteksi. Simulasi berjalan di CPU (lambat).")

print("2. Memuat model pretrained CHGNet ke " + device)
chgnet = CHGNet.load(use_device=device)

# --- SOLUSI UTAMA: Membuat kalkulator manual dengan stress=False ---
# Ini mencegah PyTorch melakukan autograd berat pada strain graph
calc = CHGNetCalculator(model=chgnet, stress=False)
atoms.set_calculator(calc)
# ------------------------------------------------------------------

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
print("   Konversi Berhasil!")
