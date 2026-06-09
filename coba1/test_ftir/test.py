import torch

# --- MODIFIKASI DISINI UNTUK GPU ---
# Mengecek apakah GPU NVIDIA (CUDA) tersedia di komputer Anda
if torch.cuda.is_available():
    device = "cuda"
    print("-> GPU terdeteksi! Simulasi akan berjalan di GPU.")
else:
    device = "cpu"
    print("cpu")
