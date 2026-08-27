#!/usr/bin/env python3
"""select_runs.py -- choose liver/lung SRA runs per BioProject for one species.

stdin : lines  Run<TAB>size_MB<TAB>BioProject<TAB>tissue   (tissue in {liver,lung})
stdout: lines  projrank<TAB>BioProject<TAB>Run<TAB>tissue
        projects ordered by total size (largest first); within a project up to 4
        runs chosen by size with the tissue rule:
          both tissues: nL>=nU -> 3 liver + 1 lung ; else 2 lung + 2 liver
          only liver -> 4 liver ; only lung -> 4 lung ; (capped to availability)
The caller tries projects in order until the primary copy (C9A) shows expression.
"""
import sys
from collections import defaultdict

def pick(runs):
    L = sorted([r for r in runs if r[2] == "liver"], key=lambda r: -r[1])
    U = sorted([r for r in runs if r[2] == "lung"],  key=lambda r: -r[1])
    if L and U:
        if len(L) >= len(U): tL, tU = 3, 1
        else: tL, tU = 2, 2
    elif L: tL, tU = 4, 0
    else:   tL, tU = 0, 4
    return L[:tL] + U[:tU]

def main():
    by_proj = defaultdict(list)
    for ln in sys.stdin:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4: continue
        run, size, proj, tis = f[0], f[1], f[2] or "NA", f[3]
        try: size = int(float(size))
        except ValueError: size = 0
        by_proj[proj].append((run, size, tis))
    # dedupe runs that appear under both liver & lung (keep one)
    order = sorted(by_proj, key=lambda p: -sum(r[1] for r in by_proj[p]))
    for rank, proj in enumerate(order, 1):
        seen = set(); uniq = []
        for r in by_proj[proj]:
            if r[0] in seen: continue
            seen.add(r[0]); uniq.append(r)
        for run, size, tis in pick(uniq):
            print(f"{rank}\t{proj}\t{run}\t{tis}")

if __name__ == "__main__":
    main()
