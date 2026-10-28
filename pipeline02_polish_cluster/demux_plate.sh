#!/usr/bin/env bash
# Demultiplex one plate, trim primers, and screen out concatemers.
#
# Usage: bash demux_plate.sh <plate> <reads.fastq.gz> [fasta_dir] [out_root] [threads] [max_len] [min_len]
#   e.g. bash demux_plate.sh 1 concat_q18.fastq.gz fasta demux 12 1500 1
# Every step uses -m so no zero-length records reach downstream tools; min_len (default 1)
# applies to the final output and can be raised to drop short junk inserts.
#
# Expects <fasta_dir>/plate_<plate>.fasta from make_plate_fastas.py.
# Writes <out_root>/plate_<plate>/final/<plate>_<row>_<sample>.fastq.gz plus counts.tsv.
set -euo pipefail

PLATE=$1; READS=$2; FASTA_DIR=${3:-fasta}; OUT_ROOT=${4:-demux}; THREADS=${5:-12}; MAXLEN=${6:-1500}; MINLEN=${7:-1}
PF=$8     # forward primer (5' -> 3' direction)
PR=`echo $9 | tr 'ACGTacgtRkYmSwHbVdDn' 'TGCAtgcyRyKmWsBdVhDn' | rev`  # reverse primer (5' -> 3' direction)

ANCHOR=8                           # must match make_plate_fastas.py --anchor
PF_REST=${PF:ANCHOR}               # primer left on the read after demux
PR_REST=${PR:0:${#PR}-ANCHOR}

OUT=$OUT_ROOT
mkdir -p "$OUT"/demux "$OUT"/trim "$OUT"/final

# Step 2: demux. The 8-nt primer anchor in each adapter keeps barcodes from matching
# inside the amplicon. (If the FASTA was built with --decoys, their output is deleted.)

echo "Demultiplexing plate $PLATE..."

cutadapt --cores="$THREADS" -m 1 -e 2 --overlap $((16 + ANCHOR)) \
  -g "file:$FASTA_DIR/plate_$PLATE.fasta" --revcomp --discard-untrimmed \
  -o "$OUT/demux/{name}.fq.gz" "$READS" > "$OUT/demux.log"
rm -f "$OUT/demux/decoy.fq.gz"

echo "Demultiplexing for plate $PLATE complete. Trimming primers and filtering concatemers..."

count() { if [[ -s $1 ]]; then echo $(( $(gzip -dc "$1" | wc -l) / 4 )); else echo 0; fi; }
printf "sample\tdemuxed\tprimers_trimmed\tfinal\n" > "$OUT/counts.tsv"

for f in "$OUT"/demux/*.fq.gz; do
  s=$(basename "$f" .fq.gz)
  # Step 3: remove the rest of each primer, anchored to the read ends.
  cutadapt --cores="$THREADS" -m 1 -e 0.15 -g "^${PF_REST}...${PR_REST}\$" \
    --discard-untrimmed -o "$OUT/trim/$s.fq.gz" "$f" > "$OUT/trim/$s.log"
  # Step 4: drop reads that still contain a primer (concatemer junction) or are too long.
  cutadapt --cores="$THREADS" -e 0.15 --overlap 18 -b "$PF" -b "$PR" --revcomp \
    --discard-trimmed -m "$MINLEN" -M "$MAXLEN" -o "$OUT/final/$s.fastq.gz" "$OUT/trim/$s.fq.gz" > "$OUT/final/$s.log"
  printf "%s\t%s\t%s\t%s\n" "$s" "$(count "$f")" "$(count "$OUT/trim/$s.fq.gz")" "$(count "$OUT/final/$s.fastq.gz")" >> "$OUT/counts.tsv"
  echo "Trimming of primers for sample_$s complete."
done

echo "plate $PLATE: $(grep 'Total reads processed' "$OUT/demux.log")"
awk -F'\t' 'NR>1{d+=$2;t+=$3;f+=$4} END{printf "  demuxed %d, primers trimmed %d, final %d reads across %d samples\n", d,t,f,NR-1}' "$OUT/counts.tsv"