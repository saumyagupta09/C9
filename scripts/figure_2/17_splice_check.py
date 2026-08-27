#!/usr/bin/env python3
"""17_splice_check.py -- compare RNA-seq splice junctions to our exon boundaries.

REPORT ONLY -- never edits exon_seqs_final. For one species, reads STAR's
aggregated splice junctions (SJ.out.tab) and our frame-exact exon coordinates
(parsed/exon_coords.tsv), and for every expected intron (gap between adjacent
exons of a copy) reports whether RNA-seq confirms it, shows a shifted boundary,
or is unsupported; also flags novel (read-supported) junctions inside the locus
not matching any expected intron (possible alternative splicing).

Appends to parsed/rnaseq_splice_check.tsv.
Usage:  python3 17_splice_check.py <RUN_DIR> <species> <SJ.tab>
"""
import sys, os
from collections import defaultdict

TOL = 5          # bp tolerance for a "shifted" boundary
MIN_NOVEL = 3    # min unique reads to flag a novel junction

def main(run, sp, sjf):
    P = os.path.join(run, "parsed")
    # our exon coords per locus  (chrom, start0, end, strand) from step 09
    ex = defaultdict(list)
    with open(os.path.join(P, "exon_coords.tsv")) as fh:
        next(fh, None)
        for ln in fh:
            lid, idx, chrom, s0, e, strand = ln.rstrip("\n").split("\t")
            ex[lid].append((int(idx), chrom, int(s0), int(e), strand))
    # locus -> (species, synteny name)
    name = {}
    with open(os.path.join(P, "master_copies.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            name[f[ci["locus_id"]]] = (f[ci["species"]], f[ci["synteny_name"]])
    # STAR junctions: chrom -> list of (istart, iend, uniq_reads)
    sj = defaultdict(list)
    if os.path.exists(sjf):
        for ln in open(sjf):
            p = ln.split("\t")
            if len(p) < 7: continue
            sj[p[0]].append((int(p[1]), int(p[2]), int(p[6])))

    rows = []
    for lid, exons in ex.items():
        if name.get(lid, ("", ""))[0] != sp: continue
        cname = name[lid][1]
        chrom = exons[0][1]
        g = sorted(((s0, e) for _, _, s0, e, _ in exons))   # by genomic position
        introns = [(g[i][1] + 1, g[i + 1][0]) for i in range(len(g) - 1)]  # (donor, acceptor) 1-based
        cstart, cend = g[0][0], g[-1][1]
        used = set()
        for k, (istart, iend) in enumerate(introns, 1):
            best = None
            for j, (sjs, sje, ur) in enumerate(sj.get(chrom, [])):
                if abs(sjs - istart) <= TOL and abs(sje - iend) <= TOL:
                    if best is None or ur > best[2]: best = (sjs, sje, ur, j)
            if best is None:
                rows.append([sp, cname, str(k), chrom, str(istart), str(iend), "-", "-", "0", "no_junction"])
            else:
                used.add(best[3])
                exact = (best[0] == istart and best[1] == iend)
                st = "confirmed" if exact and best[2] > 0 else ("shifted" if not exact else "unsupported")
                rows.append([sp, cname, str(k), chrom, str(istart), str(iend),
                             str(best[0]), str(best[1]), str(best[2]), st])
        # novel junctions inside the locus not matching expected introns
        for j, (sjs, sje, ur) in enumerate(sj.get(chrom, [])):
            if j in used or ur < MIN_NOVEL: continue
            if cstart <= sjs <= cend and cstart <= sje <= cend:
                rows.append([sp, cname, "novel", chrom, "-", "-", str(sjs), str(sje), str(ur), "novel_junction"])

    out = os.path.join(P, "rnaseq_splice_check.tsv")
    new = not os.path.exists(out)
    with open(out, "a") as fh:
        if new: fh.write("species\tcopy\tintron\tchrom\texp_donor\texp_acceptor\t"
                         "rnaseq_donor\trnaseq_acceptor\tuniq_reads\tstatus\n")
        for r in rows: fh.write("\t".join(r) + "\n")
    conf = sum(1 for r in rows if r[9] == "confirmed")
    print(f"[17] {sp}: {len(rows)} junction checks ({conf} confirmed) -> rnaseq_splice_check.tsv")

if __name__ == "__main__":
    if len(sys.argv) != 4: sys.exit("usage: 17_splice_check.py <RUN_DIR> <species> <SJ.tab>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
