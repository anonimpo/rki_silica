import torch
import ase.io
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
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

# --- KONFIGURASI CPU ---
device = "cpu"
print("-> Menjalankan simulasi di CPU.")
torch.set_num_threads(8)
print("-> Menggunakan 8 CPU threads (optimal untuk performance).")

print("2. Memuat model pretrained CHGNet ke " + device)
chgnet = CHGNet.load(use_device=device)
# ------------------------

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
