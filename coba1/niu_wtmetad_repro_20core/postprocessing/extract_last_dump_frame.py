#!/usr/bin/env python3
from pathlib import Path
import argparse

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('dump')
    ap.add_argument('--out', required=True)
    args=ap.parse_args()
    lines=Path(args.dump).read_text(errors='ignore').splitlines()
    starts=[i for i,l in enumerate(lines) if l.startswith('ITEM: TIMESTEP')]
    if not starts: raise SystemExit('No frames found')
    start=starts[-1]
    end=starts[-1+1] if False else len(lines)
    Path(args.out).write_text('\n'.join(lines[start:end])+'\n')
    print(args.out)
if __name__=='__main__': main()
