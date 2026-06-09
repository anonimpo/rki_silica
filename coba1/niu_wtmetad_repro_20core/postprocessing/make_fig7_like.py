#!/usr/bin/env python3
"""Make a Niu Fig. 7-like six-panel image from LAMMPS dump frames.
This is a visualization helper. Crystal-like Si selection is an approximate Si-Si diamond-neighbor score,
not the exact local entropy fingerprint unless MULTICOLVAR data are provided separately.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

NN_DIST=7.15*math.sqrt(3)/4  # beta-cristobalite Si-Si nearest neighbor in idealized diamond Si framework

def read_dump(path):
    lines=Path(path).read_text(errors='ignore').splitlines(); frames=[]; i=0
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            ts=int(lines[i+1].strip()); n=int(lines[i+3].strip())
            bounds=[]; j=i+5
            for _ in range(3):
                lo,hi,*_=lines[j].split(); bounds.append((float(lo),float(hi))); j+=1
            header=lines[j].split()[2:]; cols={h:k for k,h in enumerate(header)}; j+=1
            ids=[]; types=[]; xyz=[]
            for k in range(n):
                toks=lines[j+k].split()
                ids.append(int(float(toks[cols['id']])))
                types.append(int(float(toks[cols['type']])))
                xyz.append([float(toks[cols.get('x',cols.get('xu',0))]),float(toks[cols.get('y',cols.get('yu',0))]),float(toks[cols.get('z',cols.get('zu',0))])])
            frames.append((ts,np.array(ids),np.array(types),np.array(xyz),np.array([b[1]-b[0] for b in bounds])))
            i=j+n
        else:
            i+=1
    return frames

def crystal_score(si_pos, box, tol=0.35):
    n=len(si_pos); score=np.zeros(n)
    for i in range(n):
        dv=si_pos-si_pos[i]
        dv-=box*np.round(dv/box)
        d=np.linalg.norm(dv,axis=1)
        count=np.count_nonzero((d>NN_DIST-tol)&(d<NN_DIST+tol))
        # Ideal diamond Si has 4 nearest Si neighbors; map more ordered -> lower entropy-like value.
        score[i]= -0.97 - 0.78*min(count,4)/4.0
    return score

def draw_box(ax, L):
    corners=np.array([[0,0,0],[L[0],0,0],[0,L[1],0],[0,0,L[2]],[L[0],L[1],0],[L[0],0,L[2]],[0,L[1],L[2]],[L[0],L[1],L[2]]])
    edges=[(0,1),(0,2),(0,3),(1,4),(1,5),(2,4),(2,6),(3,5),(3,6),(4,7),(5,7),(6,7)]
    for a,b in edges:
        ax.plot(*zip(corners[a],corners[b]), lw=0.6, color='black', alpha=0.55)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('dump')
    ap.add_argument('--out', default='outputs/fig7_like.png')
    ap.add_argument('--every-ps', type=float, default=40.0)
    ap.add_argument('--dt-fs', type=float, default=0.25)
    ap.add_argument('--dump-stride', type=int, default=160000)
    ap.add_argument('--threshold', type=float, default=-1.16, help='show Si atoms with entropy-like score <= threshold')
    args=ap.parse_args()
    frames=read_dump(args.dump)
    if not frames: raise SystemExit('No frames found')
    # Take first 6 frames if dump already every 40 ps; otherwise spread across trajectory.
    if len(frames)>=6:
        idx=np.linspace(0,len(frames)-1,6,dtype=int)
        selected=[frames[i] for i in idx]
    else:
        selected=frames
    fig=plt.figure(figsize=(12,8))
    labels=list('ABCDEF')
    last_sc=None
    for p,fr in enumerate(selected):
        ts,ids,types,xyz,box=fr
        ax=fig.add_subplot(2,3,p+1,projection='3d')
        si=xyz[types==1]
        score=crystal_score(si,box)
        mask=score<=args.threshold
        shown=si[mask]; vals=score[mask]
        if len(shown)==0:
            shown=si; vals=score
        sc=ax.scatter(shown[:,0],shown[:,1],shown[:,2],c=vals,s=16,vmin=-1.75,vmax=-0.97,alpha=0.8)
        last_sc=sc
        draw_box(ax,box)
        ax.set_title(labels[p],loc='left',fontweight='bold',fontsize=18)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_xlim(0,box[0]); ax.set_ylim(0,box[1]); ax.set_zlim(0,box[2])
        ax.view_init(elev=17,azim=-62)
    if last_sc:
        cbar=fig.colorbar(last_sc, ax=fig.axes, shrink=0.55, pad=0.02)
        cbar.set_label('entropy-like order score')
    fig.suptitle('Fig. 7-like snapshots: only crystal-like Si shown; O/liquid-like atoms omitted',y=0.98)
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    plt.savefig(args.out,dpi=220,bbox_inches='tight')
    print(args.out)
if __name__=='__main__': main()
