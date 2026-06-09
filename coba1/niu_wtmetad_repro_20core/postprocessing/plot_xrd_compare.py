#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--before', required=True)
    ap.add_argument('--after', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--title', default='XRD before/after WTMETAD check')
    args=ap.parse_args()
    b=np.loadtxt(args.before, comments='#')
    a=np.loadtxt(args.after, comments='#')
    plt.figure(figsize=(9,5))
    plt.plot(b[:,0], b[:,1], label='before / amorphous-like')
    plt.plot(a[:,0], a[:,1]+120, label='after / beta-cristobalite-like target')
    plt.xlabel('2θ (degree), Cu Kα λ=1.5406 Å')
    plt.ylabel('normalized intensity + offset')
    plt.title(args.title)
    plt.legend()
    plt.tight_layout()
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    plt.savefig(args.out, dpi=220)
    print(args.out)
if __name__=='__main__': main()
