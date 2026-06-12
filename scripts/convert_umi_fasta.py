import sys

def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f_in, open(outfile, "w") as f_out:
    for line in f_in:
        parts = line.strip().split("\t")
        sample = parts[0]
        fwd = parts[1]
        rev = parts[2]

        rev_rc = revcomp(rev)
        full_seq = fwd + "..." + rev_rc

        f_out.write(f">{sample}\n{full_seq}\n")
