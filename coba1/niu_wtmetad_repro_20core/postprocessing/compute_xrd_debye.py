#!/usr/bin/env python3
"""Compute powder XRD from LAMMPS data/XYZ/dump using Debye scattering equation.
This is postprocessing only; it does not create MD bias.
Default uses constant scattering weights: Si=14, O=8, H=1.
"""
from __future__ import annotations
import argparse, math, re
from pathlib import Path
import numpy as np

Z = {'Si':14.0,'O':8.0,'H':1.0,1:14.0,2:8.0,3:1.0}
SYM_BY_TYPE = {1:'Si',2:'O',3:'H'}

def read_lammps_data(path: Path):
    lines = path.read_text().splitlines()
    natoms = None; box=[]; atoms_start=None
    for idx,line in enumerate(lines):
        s=line.strip()
        if s.endswith('atoms') and natoms is None:
            natoms=int(s.split()[0])
        if 'xlo xhi' in s or 'ylo yhi' in s or 'zlo zhi' in s:
            parts=s.split(); box.append(float(parts[1])-float(parts[0]))
        if s.startswith('Atoms'):
            atoms_start=idx+2; break
    if atoms_start is None: raise ValueError('Atoms section not found')
    pos=[]; types=[]
    for line in lines[atoms_start:]:
        if not line.strip() or re.match(r'^[A-Za-z]', line.strip()):
            if pos: break
            continue
        toks=line.split('#')[0].split()
        if len(toks) < 6: continue
        # atom_style charge: id type q x y z
        types.append(int(toks[1])); pos.append([float(toks[3]),float(toks[4]),float(toks[5])])
        if natoms and len(pos) >= natoms: break
    return np.array(pos,float), np.array(types,int), np.array(box,float)

def read_xyz(path: Path):
    with path.open() as f:
        n=int(f.readline().strip()); comment=f.readline().strip()
        lat=None
        m=re.search(r'Lattice="([^"]+)"', comment)
        if m:
            vals=[float(x) for x in m.group(1).split()]
            lat=np.array([vals[0],vals[4],vals[8]],float)
        syms=[]; pos=[]
        for _ in range(n):
            toks=f.readline().split(); syms.append(toks[0]); pos.append([float(toks[1]),float(toks[2]),float(toks[3])])
    types=np.array([{'Si':1,'O':2,'H':3}.get(s,0) for s in syms],int)
    if lat is None:
        arr=np.array(pos,float); lat=arr.max(axis=0)-arr.min(axis=0)
    return np.array(pos,float), types, lat

def read_lammps_dump_last(path: Path):
    text=path.read_text(errors='ignore').splitlines()
    frames=[]; i=0
    while i < len(text):
        if text[i].startswith('ITEM: TIMESTEP'):
            ts=int(text[i+1].strip()); n=int(text[i+3].strip())
            box=[]; j=i+5
            for _ in range(3):
                lo,hi,*_=text[j].split(); box.append(float(hi)-float(lo)); j+=1
            header=text[j].split()[2:]; j+=1
            cols={name:k for k,name in enumerate(header)}
            pos=[]; types=[]
            for k in range(n):
                toks=text[j+k].split()
                types.append(int(float(toks[cols.get('type',1)])))
                # prefer unwrapped/scaled? default x y z
                pos.append([float(toks[cols.get('x',cols.get('xu',2))]), float(toks[cols.get('y',cols.get('yu',3))]), float(toks[cols.get('z',cols.get('zu',4))])])
            frames.append((ts,np.array(pos),np.array(types),np.array(box)))
            i=j+n
        else:
            i+=1
    if not frames: raise ValueError('No LAMMPS dump frames found')
    return frames[-1][1], frames[-1][2], frames[-1][3]

def read_structure(path: Path):
    if path.suffix.lower()=='.xyz': return read_xyz(path)
    if 'dump' in path.name.lower() or path.suffix.lower()=='.lammpstrj': return read_lammps_dump_last(path)
    return read_lammps_data(path)

def debye(pos, types, box, two_theta, wavelength=1.5406, rcut=None, block=256, species_filter=None):
    if species_filter:
        keep=np.array([SYM_BY_TYPE.get(int(t),'X') in species_filter for t in types])
        pos=pos[keep]; types=types[keep]
    n=len(pos)
    weights=np.array([Z.get(int(t),1.0) for t in types],float)
    q=4*math.pi*np.sin(np.deg2rad(two_theta/2.0))/wavelength
    inten=np.zeros_like(q)
    if rcut is None:
        rcut=float(np.min(box)/2.0)
    for start in range(0,n,block):
        p=pos[start:start+block]
        dv=p[:,None,:]-pos[None,:,:]
        dv-=box*np.round(dv/box)
        r=np.linalg.norm(dv,axis=2).reshape(-1)
        ww=(weights[start:start+block,None]*weights[None,:]).reshape(-1)
        mask=r <= rcut
        r=r[mask]; ww=ww[mask]
        rr=np.where(r<1e-12,1e-12,r)
        window=np.sin(math.pi*rr/rcut)/(math.pi*rr/rcut)
        window[r<1e-10]=1.0
        qr=q[:,None]*rr[None,:]
        sinc=np.sin(qr)/qr
        sinc[:,r<1e-10]=1.0
        inten += np.sum(ww[None,:]*sinc*window[None,:], axis=1)
    inten=inten/n
    inten-=inten.min()
    if inten.max()>0: inten=inten/inten.max()*100.0
    return inten

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('structure')
    ap.add_argument('--out', required=True)
    ap.add_argument('--min2theta', type=float, default=10)
    ap.add_argument('--max2theta', type=float, default=90)
    ap.add_argument('--npts', type=int, default=2400)
    ap.add_argument('--lambda_', type=float, default=1.5406, dest='wavelength')
    ap.add_argument('--rcut', type=float, default=None)
    ap.add_argument('--species', default='', help='comma list e.g. Si or Si,O; blank=all')
    args=ap.parse_args()
    pos,types,box=read_structure(Path(args.structure))
    tt=np.linspace(args.min2theta,args.max2theta,args.npts)
    filt=set(args.species.split(',')) if args.species else None
    inten=debye(pos,types,box,tt,args.wavelength,args.rcut,species_filter=filt)
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    np.savetxt(args.out, np.c_[tt,inten], header='2theta_deg intensity_norm')
    print(f'Wrote {args.out} from {Path(args.structure).name}; atoms={len(pos)}, species={args.species or "all"}, box={box}')

if __name__=='__main__': main()
