#!/usr/bin/env python3
"""Build one cutadapt linked-adapter FASTA per plate from a pooled barcode sheet.

Usage:
  python make_plate_fastas.py plates_umi_demux_Run7.csv -o fasta_dir [--anchor 8] [--decoys]

Sheet: tab- or comma-delimited, no header, 4 columns:
  plate <tab> sample <tab> UDI_forward <tab> UDI_reverse   (barcodes 5'->3' as ordered)

Writes fasta_dir/plate_<P>.fasta for every plate in column 1. Each file holds:
  that plate's samples, named plate<P>_<row>_<sample> (row = position within the
  plate, so samples run more than once stay in separate files).
Optional --decoys adds every OTHER plate's combination as "decoy". In testing this
changed ~1 in 5,500 assignments, cost ~8-10x runtime, and discarded valid second
units from cross-plate concatemers, so it is off by default.

Adapter = F_barcode + first N nt of the forward primer ... last N nt of the reverse
primer (as it reads on the + strand) + R_barcode. The primer anchor stops
barcodes from matching inside the amplicon. --anchor 0 gives barcode-only adapters.
Both F and R primers should be provided from 5' to 3' direction.
"""
import argparse, collections, csv, itertools, os, sys
from Bio.Seq import Seq

# defines a function to reverse complement the primer sequence
# below function is commented out because I am trying out Biopython's Seq object for reverse complementing
#def rc(s):
#    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def edit(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sheet")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--F", required=True, help="forward primer, 5'->3'")
    ap.add_argument("--R", required=True, help="reverse primer, 5'->3'")
    ap.add_argument("--anchor", type=int, default=8, help="primer nt added next to each barcode (default 8)")
    ap.add_argument("--decoys", action="store_true", help="add other plates' combinations as 'decoy' (slow; rarely needed)")
    a = ap.parse_args()

    PRIMER_F = a.F.upper()
    PRIMER_R = a.R.upper()

    with open(a.sheet, newline="") as fh:
        text = fh.read()
    delim = "\t" if "\t" in text.splitlines()[0] else ","
    rows = [r for r in csv.reader(text.splitlines(), delimiter=delim) if r and any(x.strip() for x in r)]

    plates = collections.OrderedDict()
    combos = {}
    for n, r in enumerate(rows, 1):
        if len(r) < 4:
            sys.exit(f"ERROR line {n}: expected 4 columns, got {len(r)}")
        plate, sample, fwd, rev = (x.strip() for x in r[:4])
        fwd, rev = fwd.upper(), rev.upper()
        for bc in (fwd, rev):
            if set(bc) - set("ACGT"):
                sys.exit(f"ERROR line {n}: non-ACGT barcode {bc}")
        if (fwd, rev) in combos:
            sys.exit(f"ERROR line {n}: combination reused (also {combos[(fwd, rev)]})")
        plates.setdefault(plate, [])
        label = f"plate{plate}_{len(plates[plate]) + 1:02d}_{sample}"
        combos[(fwd, rev)] = label
        plates[plate].append((label, sample, fwd, rev))

    fa, ra = PRIMER_F[: a.anchor], Seq(PRIMER_R).reverse_complement()[len(PRIMER_R) - a.anchor:]
    adapter = lambda f, r: f"{f}{fa}...{ra}{Seq(r).reverse_complement()}"
    #adapter = lambda f, r: f"{f}{fa}...{ra}{rc(r)}" #commented out to use biopython reverse complement function

    # Report
    allF = sorted({f for f, _ in combos}); allR = sorted({r for _, r in combos})
    print(f"{len(rows)} combinations, {len(plates)} plates, {len(allF)} F / {len(allR)} R barcodes", file=sys.stderr)
    print(f"min pairwise edit distance: F {min(edit(x, y) for x, y in itertools.combinations(allF, 2))}, "
          f"R {min(edit(x, y) for x, y in itertools.combinations(allR, 2))}", file=sys.stderr)
    for p, v in plates.items():
        dup = [s for s, c in collections.Counter(x[1] for x in v).items() if c > 1]
        note = f"; replicated samples (separate files by row): {', '.join(dup)}" if dup else ""
        print(f"plate {p}: {len(v)} samples{note}", file=sys.stderr)
    shared = collections.defaultdict(set)
    for p, v in plates.items():
        for x in v:
            shared[x[1]].add(p)
    multi = {s: sorted(ps) for s, ps in shared.items() if len(ps) > 1}
    if multi:
        print("samples replicated across plates: " + ", ".join(f"{s} ({'/'.join(ps)})" for s, ps in multi.items()), file=sys.stderr)

    os.makedirs(a.outdir, exist_ok=True)
    for p, v in plates.items():
        path = os.path.join(a.outdir, f"plate_{p}.fasta")
        with open(path, "w") as out:
            for label, _, f, r in v:
                out.write(f">{label}\n{adapter(f, r)}\n")
            if a.decoys:
                for q, w in plates.items():
                    if q != p:
                        for _, _, f, r in w:
                            out.write(f">decoy\n{adapter(f, r)}\n")
        print(f"wrote {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
