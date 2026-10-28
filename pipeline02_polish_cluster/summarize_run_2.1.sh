#!/usr/bin/env bash
#
# summarize_run.sh
#
# Summarize an amplicon sequencing run from a long-form OTU table plus a BLAST6
# table. Taxonomy is derived internally (no pre-built accession/lineage files
# needed): BLAST hits -> accessions -> taxids (acc2taxid) -> taxonkit lineage +
# LCA. Optionally collapses OTUVs into taxonomic OTUs.
#
# ---------------------------------------------------------------------------
# TERMS
#   OTUV  Operational Taxonomic Unit Variant = a haplotype / variant sequence
#         (the cluster IDs in the input, e.g. cluster_15;size=4). OTUVs reflect
#         fine-scale genetic diversity and are the atomic unit of the counts.
#   OTU   Operational Taxonomic Unit = OTUVs merged by a criterion. Here the
#         criterion is shared best-hit SPECIES (from BLAST). OTUVs with no
#         resolved species do not merge; each passes through as its own OTU
#         labelled "Unassigned (<OTUV>)".
#   NTC   No-template control sample(s): contain no real template but still
#         yield reads (index hopping, contamination, carryover). Used to
#         estimate and subtract background.
#   LCA   Lowest Common Ancestor of the taxids of all coverage-passing BLAST
#         hits for an OTUV, reported at its finest resolved rank.
# ---------------------------------------------------------------------------
# PIPELINE ORDER (fixed)
#   1. NTC subtraction   (--ntc): per OTUV, take the MAX read count across the
#      named NTC samples and subtract it from that OTUV in every sample.
#      Negative results clamp to 0.
#   2. Low-count filter  (--lowcount): per (sample, OTUV), if the NTC-adjusted
#      read count is below the threshold, set it to 0.
#   3. mindepth check    (--mindepth): flag samples whose total reads AFTER
#      NTC+lowcount but BEFORE taxonomy fall below the threshold. Taxonomy loss
#      reflects primer non-specificity, so it is excluded from this check.
#   4. Taxonomy merge    (--taxonomy): optionally collapse OTUVs into OTUs by
#      shared best-hit species, optionally restricted to phyla (--phylum).
#   5. Reports.
# ---------------------------------------------------------------------------
# INPUTS
#   -c FILE   OTU table, LONG form: Sample <tab> OTUV <tab> Reads  (header row)
#   -b FILE   BLAST6 table. Expected custom outfmt (matches:
#               '6 qseqid sseqid pident length qlen slen qstart qend sstart
#                  send evalue bitscore'):
#               1 qseqid(OTUV) 2 sseqid 3 pident 4 aln_len 5 qlen 6 slen
#               ... 12 bitscore
#             Query coverage = aln_len/qlen (col4/col5). Best hit = max bitscore.
#
# TAXONOMY TOOL PATHS (no defaults assumed)
#   --acc2taxid FILE   NCBI accession2taxid map. Matched on accession.version
#                      in col2 -> taxid in col3 (header row auto-skipped);
#                      falls back to col1->col2.
#   --taxonkit PATH    path to the taxonkit binary
#   --taxdump DIR      path to the NCBI taxdump directory
#
# FILTERING / OPTIONS
#   --ntc=LIST         comma list of NTC samples (per-OTUV max subtraction)
#   --lowcount=N       per (sample,OTUV) reads < N -> 0            (default 0/off)
#   --mindepth=N       flag samples below N reads post-NTC+lowcount, pre-tax
#                      (default 1000)
#   --mincov=F         min query coverage (0-1) for a hit to be eligible as
#                      best hit and to feed the LCA               (default 0.25)
#   --taxonomy         collapse OTUVs into OTUs by shared best-hit species
#   --phylum=LIST      Phylum can be used with and without a list. When used without
#                      a list, then keep only OTUVs that have a defined phylum in
#                      their hits. Or, specify a comma-separated list - and keep 
#                      only OTUVs whose best-hit phylum is in this comma list
#                      (case-insensitive). OTUVs with no/blank/unresolved phylum 
#                      are dropped too. Applies to ALL outputs and works with 
#                      or without --taxonomy.
#   -o PREFIX          output filename prefix
#   --quiet            suppress progress messages (otherwise printed to stderr)
#
# OUTPUTS (prefixed with -o PREFIX)
#   {p}samples_report.txt   Per sample: reads after filtering + top OTU
#                           (species, with min-max %ID/%cov) or top OTUV.
#                           Retention + NTC footer + mindepth summary.
#   {p}otuvs_per_sample.txt Per sample: OTUV (and OTU, if --taxonomy) counts,
#                           surviving list, and read-loss split into
#                           NTC / lowcount / taxonomy (reads and % of raw).
#   {p}otu_summary.txt      Detailed OTU/OTUV report: per merged OTU, the member
#                           OTUVs with individual %ID/%cov/LCA/confidence, plus
#                           aggregate ranges; and contaminant screening.
#   {p}otu_counts.tsv       (--taxonomy only) wide matrix samples x OTU(species),
#                           filtered cells NA, NTC samples excluded.
#   {p}otuv_counts.tsv      wide matrix samples x OTUV, filtered cells NA,
#                           NTC excluded.
#   {p}otuv_long.tsv        LONG: Sample OTUV Reads pident qcov phylum class
#                           order family subfamily genus species LCA n_hits
#                           n_ties confidence
#
# CONFIDENCE TAG (per OTUV; in otuv_long.tsv):
#   concordant       LCA resolves as fine as the best hit
#   ambiguous-<rank> best hit reaches species but LCA backs off to <rank>
#   LCA-divergent    LCA sits on a different lineage than the best hit
#   LCA-unresolved   LCA did not resolve
#   no-hit           no BLAST hit passed coverage
#
set -euo pipefail

COUNTS=""; B6=""
ACC2TAXID=""; TAXONKIT=""; TAXDUMP=""
NTC_LIST=""; LOWCOUNT=0; MINDEPTH=1000; MINCOV=0.25
DO_TAX=0; PHYLUM_LIST=""; PHYLUM_SET=0
PREFIX=""; QUIET=0

usage() { sed -n '2,110p' "$0"; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
    case "$1" in
        -c) COUNTS="$2"; shift 2 ;;
        -b) B6="$2"; shift 2 ;;
        --acc2taxid=*) ACC2TAXID="${1#--acc2taxid=}"; shift ;;
        --acc2taxid)   ACC2TAXID="$2"; shift 2 ;;
        --taxonkit=*)  TAXONKIT="${1#--taxonkit=}"; shift ;;
        --taxonkit)    TAXONKIT="$2"; shift 2 ;;
        --taxdump=*)   TAXDUMP="${1#--taxdump=}"; shift ;;
        --taxdump)     TAXDUMP="$2"; shift 2 ;;
        --ntc=*)       NTC_LIST="${1#--ntc=}"; shift ;;
        --ntc)         NTC_LIST="$2"; shift 2 ;;
        --lowcount=*)  LOWCOUNT="${1#--lowcount=}"; shift ;;
        --lowcount)    LOWCOUNT="$2"; shift 2 ;;
        --mindepth=*)  MINDEPTH="${1#--mindepth=}"; shift ;;
        --mindepth)    MINDEPTH="$2"; shift 2 ;;
        --mincov=*)    MINCOV="${1#--mincov=}"; shift ;;
        --mincov)      MINCOV="$2"; shift 2 ;;
        --taxonomy)    DO_TAX=1; shift ;;
        --phylum=*)    PHYLUM_SET=1; PHYLUM_LIST="${1#--phylum=}"; shift ;;
        --phylum)      PHYLUM_SET=1; PHYLUM_LIST=""; shift ;;
        -o) PREFIX="$2"; shift 2 ;;
        --quiet) QUIET=1; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

[ -z "$COUNTS" ] && { echo "Missing -c OTU table." >&2; usage 1; }
[ -f "$COUNTS" ] || { echo "File not found: $COUNTS" >&2; exit 1; }
[ -z "$B6" ] && { echo "Missing -b BLAST6 table." >&2; usage 1; }
[ -f "$B6" ] || { echo "File not found: $B6" >&2; exit 1; }
for req in ACC2TAXID TAXONKIT TAXDUMP; do
    eval "v=\$$req"
    if [ -z "$v" ]; then
        lc=$(printf '%s' "$req" | tr 'A-Z' 'a-z')
        echo "Missing --$lc path (no default)." >&2; exit 1
    fi
done
[ -f "$ACC2TAXID" ] || { echo "acc2taxid not found: $ACC2TAXID" >&2; exit 1; }
[ -x "$TAXONKIT" ] || command -v "$TAXONKIT" >/dev/null 2>&1 || { echo "taxonkit not executable: $TAXONKIT" >&2; exit 1; }
[ -d "$TAXDUMP" ] || { echo "taxdump dir not found: $TAXDUMP" >&2; exit 1; }
P="$PREFIX"

# Progress logging -> stderr (never pollutes the output files on stdout).
# Silence with --quiet, or by redirecting stderr.
log() { [ "$QUIET" = 1 ] && return 0; printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

TMPDIR="$(mktemp -d)"
trap 'rm -f "$TMPDIR"/*; rmdir "$TMPDIR" 2>/dev/null || true' EXIT

log "summarize_run: starting"
log "  counts:  $COUNTS"
log "  blast6:  $B6"
log "  options: ntc=[${NTC_LIST:-none}] lowcount=$LOWCOUNT mindepth=$MINDEPTH mincov=$MINCOV taxonomy=$DO_TAX phylum=[${PHYLUM_LIST:-$( [ "$PHYLUM_SET" = 1 ] && echo any || echo none )}]"

# ======================================================================
# STEP 1: parse BLAST6 -> per-OTUV best hit + coverage-passing hit set.
#   BESTHIT : OTUV  best_accession  pident  qcov  bitscore
#   HITACC  : OTUV  accession        (one row per coverage-passing hit)
# ======================================================================
log "Step 1/6: parsing BLAST6, selecting best hits (coverage >= $MINCOV)..."
BESTHIT="$TMPDIR/besthit.tsv"
HITACC="$TMPDIR/hitacc.tsv"
awk -v FS='\t' -v OFS='\t' -v MINCOV="$MINCOV" -v BESTHIT="$BESTHIT" -v HITACC="$HITACC" '
function acc_of(s,   n,p,i){ if(index(s,"|")>0){ n=split(s,p,"|"); for(i=1;i<=n;i++) if(p[i] ~ /^[A-Z][A-Z]?[0-9][0-9]*(\.[0-9]+)?$/) return p[i]; if(n>=4)return p[4] } return s }
{
    q=$1; acc=acc_of($2); pid=$3+0; alen=$4+0; qlen=$5+0; bits=$12+0
    cov=(qlen>0)?alen/qlen:0
    if (cov < MINCOV) next
    nhits[q]++                                    # coverage-passing hits
    k=q SUBSEP acc
    if (!(k in seen)) { seen[k]=1; print q, acc > HITACC }
    if (!(q in bestb) || bits > bestb[q]) {
        bestb[q]=bits; besta[q]=acc; bestp[q]=pid; bestc[q]=cov; ntie[q]=1
    } else if (bits == bestb[q]) {
        ntie[q]++                                 # another hit at the max bitscore
    }
}
END {
    for (q in bestb)
        printf "%s\t%s\t%.3f\t%.4f\t%s\t%d\t%d\n", q, besta[q], bestp[q], bestc[q], bestb[q], nhits[q], ntie[q] > BESTHIT
}
' "$B6"
log "  ...$(wc -l < "$BESTHIT" | tr -d ' ') OTUVs with a best hit; $(cut -f1 "$HITACC" | sort -u | wc -l | tr -d ' ') OTUVs have coverage-passing hits for LCA"

# ======================================================================
# STEP 2: accession -> taxid via acc2taxid.
# We need taxids for: (a) every best-hit accession, (b) every hit accession
# (for LCA). Build one accession list, look up once.
# acc2taxid assumed 2+ cols with accession.version in col1 (or col2) and
# taxid in the last-but-one/con configuration used previously:
#   previously matched: ($2 in a){print $2"\t"$3}  i.e. col2=acc.ver col3=taxid
# We match on accession.version against col2, emit col3 as taxid. Fallback:
# if col1 looks like the accession, use col1->col2.
# ======================================================================
ACCLIST="$TMPDIR/acclist.txt"
cut -f2 "$HITACC" | sort -u > "$ACCLIST"
log "Step 2/6: resolving $(wc -l < "$ACCLIST" | tr -d ' ') unique accessions to taxids (acc2taxid may be large; this can take a moment)..."

ACC2TAX="$TMPDIR/acc2tax.tsv"   # accession \t taxid
awk -v FS='\t' -v OFS='\t' '
NR==FNR { want[$1]=1; next }
FNR==1 && ($1=="accession" || $2=="accession.version") { next }   # skip header
{
    # try col2=acc.ver col3=taxid (NCBI style used before)
    if ($2 in want) { print $2, $3; next }
    # fallback col1=acc col2=taxid
    if ($1 in want) { print $1, $2 }
}
' "$ACCLIST" "$ACC2TAXID" | sort -u > "$ACC2TAX"
log "  ...mapped $(wc -l < "$ACC2TAX" | tr -d ' ') accessions to taxids"

# ======================================================================
# STEP 3: taxonkit.
# (a) Best-hit lineage: taxid -> phylum,class,order,family,subfamily,genus,species
# (b) LCA per OTUV: gather taxids of all coverage-passing hits, taxonkit lca,
#     then lineage|reformat to get finest-resolution name.
# ======================================================================
log "Step 3/6: running taxonkit (best-hit lineages + per-OTUV LCA)..."
# map accession->taxid into awk-friendly file, then join to besthit and hitacc
BEST_TAXID="$TMPDIR/best_taxid.tsv"     # OTUV taxid pident qcov nhits nties
awk -v FS='\t' -v OFS='\t' '
NR==FNR{ a2t[$1]=$2; next }
{ q=$1; acc=$2; pid=$3; cov=$4; nh=$6; nt=$7; t=(acc in a2t)?a2t[acc]:""; print q, t, pid, cov, nh, nt }
' "$ACC2TAX" "$BESTHIT" > "$BEST_TAXID"

# Best-hit lineage via taxonkit (only rows with a taxid)
BEST_LINEAGE="$TMPDIR/best_lineage.tsv"  # taxid P C O F SF G S
cut -f2 "$BEST_TAXID" | awk 'NF>0' | sort -u > "$TMPDIR/best_taxids.txt"
if [ -s "$TMPDIR/best_taxids.txt" ]; then
    log "  ...taxonkit lineage for $(wc -l < "$TMPDIR/best_taxids.txt" | tr -d ' ') best-hit taxids"
    "$TAXONKIT" lineage --data-dir "$TAXDUMP" "$TMPDIR/best_taxids.txt" \
      | "$TAXONKIT" reformat --data-dir "$TAXDUMP" -f "{p}\t{c}\t{o}\t{f}\t{S}\t{g}\t{s}" \
      | awk -v FS='\t' -v OFS='\t' '
          function clean(x){ if(x ~ /^unclassified /) return ""; return x }
          { print $1, clean($3),clean($4),clean($5),clean($6),clean($7),clean($8),clean($9) }' \
      > "$BEST_LINEAGE" 2> "$TMPDIR/tk_best.err" || { echo "taxonkit (best lineage) failed:" >&2; cat "$TMPDIR/tk_best.err" >&2; exit 1; }
fi

# LCA per OTUV: build "OTUV<tab>taxid,taxid,..." then taxonkit lca
HIT_TAXID="$TMPDIR/hit_taxid.tsv"
awk -v FS='\t' -v OFS='\t' '
NR==FNR{ a2t[$1]=$2; next }
{ q=$1; acc=$2; if(acc in a2t && a2t[acc]!="") print q, a2t[acc] }
' "$ACC2TAX" "$HITACC" > "$HIT_TAXID"

# group taxids per OTUV, comma-separated
OTUV_TAXIDS="$TMPDIR/otuv_taxids.tsv"   # OTUV  t1,t2,t3
awk -v FS='\t' -v OFS='\t' '
{ if($1 in g) g[$1]=g[$1]","$2; else g[$1]=$2 }
END{ for(q in g) print q, g[q] }
' "$HIT_TAXID" > "$OTUV_TAXIDS"

# taxonkit lca on the comma list (2nd field)
OTUV_LCA="$TMPDIR/otuv_lca.tsv"   # OTUV  lca_taxid
if [ -s "$OTUV_TAXIDS" ]; then
    log "  ...taxonkit lca for $(wc -l < "$OTUV_TAXIDS" | tr -d ' ') OTUVs"
    "$TAXONKIT" lca --data-dir "$TAXDUMP" -i 2 -s "," "$OTUV_TAXIDS" 2> "$TMPDIR/tk_lca.err" \
      | awk -v FS='\t' -v OFS='\t' '{ print $1, $NF }' > "$OTUV_LCA" \
      || { echo "taxonkit lca failed:" >&2; cat "$TMPDIR/tk_lca.err" >&2; exit 1; }
fi

# LCA taxid -> finest-resolution name (species if present, else genus, ... else phylum)
LCA_NAME="$TMPDIR/lca_name.tsv"   # OTUV  lca_name
if [ -s "$OTUV_LCA" ]; then
    # unique lca taxids -> reformat, then choose finest non-empty rank
    cut -f2 "$OTUV_LCA" | awk 'NF>0 && $1!=0' | sort -u > "$TMPDIR/lca_taxids.txt"
    if [ -s "$TMPDIR/lca_taxids.txt" ]; then
        "$TAXONKIT" lineage --data-dir "$TAXDUMP" "$TMPDIR/lca_taxids.txt" \
          | "$TAXONKIT" reformat --data-dir "$TAXDUMP" -f "{p}\t{c}\t{o}\t{f}\t{S}\t{g}\t{s}" \
          | awk -v FS='\t' -v OFS='\t' '
                function clean(x){ if(x ~ /^unclassified /) return ""; return x }
                {
                    taxid=$1; P=clean($3);C=clean($4);O=clean($5);F=clean($6);SF=clean($7);G=clean($8);S=clean($9)
                    name = (S!="")?S : (G!="")?G : (SF!="")?SF : (F!="")?F : (O!="")?O : (C!="")?C : (P!="")?P : ""
                    print taxid, name
                }' > "$TMPDIR/lca_taxid_name.tsv"
    fi
    awk -v FS='\t' -v OFS='\t' '
        NR==FNR{ nm[$1]=$2; next }
        { t=$2; print $1, (t in nm && nm[t]!="")?nm[t]:"Unavailable" }
    ' "$TMPDIR/lca_taxid_name.tsv" "$OTUV_LCA" > "$LCA_NAME"
fi

# ======================================================================
# STEP 4: assemble per-OTUV taxonomy annotation file consumed by reporting.
#   OTUVANNO: OTUV pident qcov P C O F SF G S LCAname best_accession
# ======================================================================
log "Step 4/6: assembling per-OTUV taxonomy annotation (best hit + LCA + confidence)..."
OTUVANNO="$TMPDIR/otuv_anno.tsv"
awk -v FS='\t' -v OFS='\t' '
# best_lineage: taxid P C O F SF G S
FILENAME==ARGV[1] { bl[$1]=$2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8; next }
# best_taxid: OTUV taxid pident qcov nhits nties
FILENAME==ARGV[2] { bt_taxid[$1]=$2; bt_pid[$1]=$3; bt_cov[$1]=$4; bt_nh[$1]=$5; bt_nt[$1]=$6; next }
# besthit: OTUV acc pident qcov bits nhits nties  (for accession)
FILENAME==ARGV[3] { bacc[$1]=$2; next }
# lca_name: OTUV name
FILENAME==ARGV[4] { lca[$1]=$2; next }
# rank depth of a lineage name by matching it against the 7 rank fields
function finest_depth(P,C,O,Fm,SF,G,S){
    if(S!="")return 7; if(G!="")return 6; if(SF!="")return 5; if(Fm!="")return 4
    if(O!="")return 3; if(C!="")return 2; if(P!="")return 1; return 0
}
# which rank does the LCA name correspond to, given the best-hit lineage?
function lca_depth(name,P,C,O,Fm,SF,G,S){
    if(name==""||name=="Unavailable")return 0
    if(name==S)return 7; if(name==G)return 6; if(name==SF)return 5; if(name==Fm)return 4
    if(name==O)return 3; if(name==C)return 2; if(name==P)return 1
    return -1   # LCA name not found in best-hit lineage (different branch)
}
END {
    for (q in bacc) {
        t=bt_taxid[q]
        lin=(t in bl)?bl[t]:(OFS OFS OFS OFS OFS OFS)
        nlin=split(lin, lf, OFS); for(i=nlin+1;i<=7;i++) lf[i]=""
        pid=(q in bt_pid)?bt_pid[q]:""
        cov=(q in bt_cov)?bt_cov[q]:""
        nh=(q in bt_nh)?bt_nh[q]:0
        nt=(q in bt_nt)?bt_nt[q]:0
        lc=(q in lca)?lca[q]:"Unavailable"

        # divergence tag: compare best-hit finest rank vs LCA finest rank
        bd=finest_depth(lf[1],lf[2],lf[3],lf[4],lf[5],lf[6],lf[7])
        ld=lca_depth(lc,lf[1],lf[2],lf[3],lf[4],lf[5],lf[6],lf[7])
        if (bd==0) tag="no-hit"
        else if (lc=="Unavailable") tag="LCA-unresolved"
        else if (ld==-1) tag="LCA-divergent"          # LCA on a different lineage than best hit
        else if (ld>=bd) tag="concordant"             # LCA as fine as (or finer than) best hit
        else {
            rn[7]="species";rn[6]="genus";rn[5]="subfamily";rn[4]="family";rn[3]="order";rn[2]="class";rn[1]="phylum"
            tag="ambiguous-" rn[ld]                   # LCA backs off to this rank
        }

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
            q, pid, cov, lf[1],lf[2],lf[3],lf[4],lf[5],lf[6],lf[7], lc, bacc[q], nh OFS nt, tag
    }
}
' "$BEST_LINEAGE" "$BEST_TAXID" "$BESTHIT" "$LCA_NAME" > "$OTUVANNO"

# ======================================================================
# STEP 5: NTC subtraction + lowcount (same engine as before), on -c long form.
# ======================================================================
log "Step 5/6: applying NTC subtraction and low-count filter..."
ADJUSTED="$TMPDIR/adjusted.tsv"
PERSAMPLE="$TMPDIR/persample.tsv"
OTUVSTATS="$TMPDIR/otuvstats.tsv"
awk -v FS='\t' -v OFS='\t' -v NTC_LIST="$NTC_LIST" -v LOWCOUNT="$LOWCOUNT" \
    -v PERSAMPLE="$PERSAMPLE" -v OTUVSTATS="$OTUVSTATS" '
BEGIN{ nn=split(NTC_LIST,a,","); for(i=1;i<=nn;i++){gsub(/^[ \t]+|[ \t]+$/,"",a[i]); if(a[i]!=""){is_ntc[a[i]]=1; ord[++m]=a[i]}} LOWCOUNT+=0 }
NR==1{ next }
{ n++; S[n]=$1;O[n]=$2;R[n]=$3+0; tot[$2]+=R[n]; seen[$2]=1
  if($1 in is_ntc){ if(!($2 in mx)||R[n]>mx[$2])mx[$2]=R[n]; nr[$2 SUBSEP $1]=R[n]; nt[$2]+=R[n] } }
END{
  print "Sample","OTUV","Reads"
  for(i=1;i<=n;i++){ r=R[i]; raw[S[i]]+=r
    sa=(O[i] in mx)?mx[O[i]]:0; rn=r-sa; if(rn<0)rn=0; pn[S[i]]+=rn; rem[S[i]]+=(r-rn)
    rl=rn; if(LOWCOUNT>0 && rl<LOWCOUNT)rl=0; pl[S[i]]+=rl; lrem[S[i]]+=(rn-rl)
    print S[i],O[i],rl }
  for(s in raw) print s,raw[s]+0,pn[s]+0,pl[s]+0,rem[s]+0,lrem[s]+0 > PERSAMPLE
  h="OTUV\ttotal_reads\tntc_reads_total\tntc_fraction"; for(k=1;k<=m;k++)h=h OFS ("ntc:" ord[k]); print h > OTUVSTATS
  for(o in seen){ tt=tot[o]+0; ntt=nt[o]+0; fr=(tt>0)?ntt/tt:0
    line=o OFS tt OFS ntt OFS sprintf("%.4f",fr)
    for(k=1;k<=m;k++){ key=o SUBSEP ord[k]; line=line OFS ((key in nr)?nr[key]:0) }
    print line > OTUVSTATS }
}
' "$COUNTS" > "$ADJUSTED"
[ -f "$PERSAMPLE" ] || : > "$PERSAMPLE"
[ -f "$OTUVSTATS" ] || : > "$OTUVSTATS"

# ======================================================================
# STEP 6: reports.
# ======================================================================
log "Step 6/6: writing reports$( [ "$DO_TAX" = 1 ] && printf ' (taxonomy-merged)'; [ "$PHYLUM_SET" = 1 ] && printf ' (phylum-filtered)' )..."
awk -v FS='\t' -v OFS='\t' -v P="$P" -v NTC_LIST="$NTC_LIST" \
    -v LOWCOUNT="$LOWCOUNT" -v MINDEPTH="$MINDEPTH" -v MINCOV="$MINCOV" \
    -v DO_TAX="$DO_TAX" -v PHYLUM_SET="$PHYLUM_SET" -v PHYLUM_LIST="$PHYLUM_LIST" \
    -v ADJ="$ADJUSTED" -v PERSAMPLE="$PERSAMPLE" -v OTUVSTATS="$OTUVSTATS" -v ANNO="$OTUVANNO" '
BEGIN{
  nn=split(NTC_LIST,a,","); for(i=1;i<=nn;i++){gsub(/^[ \t]+|[ \t]+$/,"",a[i]); if(a[i]!=""){is_ntc[a[i]]=1; ord[++m]=a[i]}}
  MINDEPTH+=0

  # annotation: OTUV pid cov P C O F SF G S LCA acc nhits nties tag
  while((getline line < ANNO)>0){ ns=split(line,f,"\t"); o=f[1]
    an_pid[o]=f[2]; an_cov[o]=f[3]; an_P[o]=f[4]; an_C[o]=f[5]; an_O[o]=f[6]; an_F[o]=f[7]; an_SF[o]=f[8]; an_G[o]=f[9]; an_S[o]=f[10]; an_LCA[o]=f[11]; an_acc[o]=f[12]; an_nh[o]=f[13]; an_nt[o]=f[14]; an_tag[o]=f[15] }
  close(ANNO)

  # per-sample stage totals
  while((getline line < PERSAMPLE)>0){ ns=split(line,f,"\t"); s=f[1]
    raw_tot[s]=f[2]+0; pn[s]=f[3]+0; pl[s]=f[4]+0; ntc_rem[s]=f[5]+0; low_rem[s]=f[6]+0; seen_s[s]=1 }
  close(PERSAMPLE)

  # Build phylum include-set (lowercased). PHYLUM_SET==1 with empty list ->
  # "any defined phylum"; PHYLUM_SET==0 -> no restriction. This is a GLOBAL
  # gate: an OTUV is dropped from ALL outputs unless its best-hit phylum is
  # named (strict: blank / no-hit / unresolved phylum are dropped too).
  n_incl=0
  if (PHYLUM_LIST!="") { np=split(PHYLUM_LIST,pp,","); for(i=1;i<=np;i++){gsub(/^[ \t]+|[ \t]+$/,"",pp[i]); if(pp[i]!=""){incl[tolower(pp[i])]=1; n_incl++}} }

  # adjusted long form (apply phylum gate here so nothing downstream sees
  # out-of-phylum OTUVs)
  getline < ADJ
  while((getline line < ADJ)>0){ ns=split(line,f,"\t"); s=f[1];o=f[2];r=f[3]+0
    seen_s[s]=1
    if(r>0){
      if (PHYLUM_SET==1) {
        phy=an_P[o]
        if (phy=="") continue                                   # strict: drop blank/no-hit
        if (n_incl>0 && !(tolower(phy) in incl)) continue        # drop other phyla
      }
      ov_r[s SUBSEP o]=r; ov_n[s]++; all_ov[o]=1; ovg_tot[s]+=r
      if(r>best_r[s]){best_r[s]=r; best_o[s]=o} } }
  close(ADJ)

  # ---- Taxonomy merge layer (only when DO_TAX==1) --------------------
  # ov_r[] is already phylum-gated above, so merging just groups survivors.
  if (DO_TAX==1) {
    for (k in ov_r) {
      split(k,kk,SUBSEP); s=kk[1]; o=kk[2]; r=ov_r[k]
      sp=an_S[o]
      phy=an_P[o]
      # OTU label: species if resolved, else Unassigned (OTUV)
      if (sp!="") otu=sp; else otu="Unassigned (" o ")"
      gk = s SUBSEP otu
      otu_r[gk]+=r
      if (!(gk in otu_seen)) { otu_seen[gk]=1; otu_n[s]++; all_otu[otu]=1; otu_phy[otu]=phy }
      otu_ovn[gk]++
      otu_tot_s[s]+=r
      if (r>otu_best_r[s]){ otu_best_r[s]=r; otu_best[s]=otu }
      otu_members[gk] = otu_members[gk] (otu_members[gk]==""?"":",") o
      if (an_pid[o]!="" && an_pid[o]!="NA") {
        v=an_pid[o]+0
        if (!(otu in otu_pidmin) || v<otu_pidmin[otu]) otu_pidmin[otu]=v
        if (!(otu in otu_pidmax) || v>otu_pidmax[otu]) otu_pidmax[otu]=v
      }
      if (an_cov[o]!="" && an_cov[o]!="NA") {
        v=100*(an_cov[o]+0)
        if (!(otu in otu_covmin) || v<otu_covmin[otu]) otu_covmin[otu]=v
        if (!(otu in otu_covmax) || v>otu_covmax[otu]) otu_covmax[otu]=v
      }
      if (!(o in otu_glob_seen)) { otu_glob_seen[o]=1
        otu_glob_members[otu] = otu_glob_members[otu] (otu_glob_members[otu]==""?"":",") o }
    }
  }

  # otuv stats (contam)
  getline < OTUVSTATS
  while((getline line < OTUVSTATS)>0){ ns=split(line,f,"\t"); o=f[1]
    ovs_tot[o]=f[2]+0; ovs_ntct[o]=f[3]+0; ovs_frac[o]=f[4]+0
    for(k=1;k<=m;k++) ovs_ntc[o SUBSEP ord[k]]=f[4+k]+0 }
  close(OTUVSTATS)

  # sample order
  ns2=0; for(s in seen_s) sa[++ns2]=s
  for(i=2;i<=ns2;i++){k=sa[i];j=i-1;while(j>=1&&sa[j]>k){sa[j+1]=sa[j];j--}sa[j+1]=k}
  # display total: OTU-level in taxonomy mode (only reads in a reported OTU),
  # else phylum-gated OTUV-level. mindepth always uses pl[] (post-NTC+lowcount,
  # pre-taxonomy AND pre-phylum: it measures sequencing depth, not assignment).
  for(i=1;i<=ns2;i++){ s=sa[i]; disp[s]=(DO_TAX==1)?(otu_tot_s[s]+0):(ovg_tot[s]+0); lowdep[s]=(pl[s]+0<MINDEPTH)?1:0 }

  # helper for "ID/cov" string
  # (inlined below)

  # ===== FILE 1: samples_report.txt =====
  f1=P "samples_report.txt"
  print "Sample sequencing report" > f1; print "=========================" > f1; print "" > f1
  if(NTC_LIST!="") print "NTC samples: " NTC_LIST > f1
  if(LOWCOUNT+0>0) print "Low-count threshold: " (LOWCOUNT+0) " reads/OTUV" > f1
  print "Min-depth threshold: " MINDEPTH " reads (post-NTC+lowcount)" > f1
  print "Min coverage for hits: " (MINCOV+0) > f1
  print "" > f1; print "Per-sample summary:" > f1; print "-------------------" > f1
  traw=0; tdisp=0; nlow=0
  for(i=1;i<=ns2;i++){ s=sa[i]; traw+=raw_tot[s]; tdisp+=disp[s]
    tg=(s in is_ntc)?"  [NTC]":""; dt=(lowdep[s]==1)?sprintf("  [LOW DEPTH: %d reads]",pl[s]+0):""
    if(lowdep[s]==1)nlow++
    if(disp[s]+0==0){ print s tg dt ": All reads from sample were filtered" > f1 }
    else if(DO_TAX==1){ otu=otu_best[s]
      rng=idcov_range(otu)
      printf "%s%s%s: %d reads | top OTU: %s (%s)\n", s, tg, dt, disp[s], otu, rng > f1 }
    else{ o=best_o[s]
      idcov=(an_pid[o]!=""&&an_pid[o]!="NA")?sprintf("%.0f%% ID, %.0f%% cov", an_pid[o]+0, 100*(an_cov[o]+0)):"no hit"
      lc=(o in an_LCA && an_LCA[o]!="")?an_LCA[o]:"Unavailable"
      printf "%s%s%s: %d reads | top OTUV: %s (%s; LCA: %s)\n", s, tg, dt, disp[s], o, idcov, lc > f1 }
  }
  print "" > f1; print "Read retention:" > f1; print "---------------" > f1
  printf "Total raw reads:      %d\n", traw > f1
  printf "Total retained reads: %d\n", tdisp > f1
  if(traw>0) printf "Overall retention:    %.1f%%\n", 100.0*tdisp/traw > f1
  printf "Samples below min-depth (%d): %d of %d\n", MINDEPTH, nlow, ns2 > f1
  print "" > f1; print "NTC subtraction (reads removed per sample):" > f1; print "-------------------------------------------" > f1
  if(NTC_LIST==""){ print "  (no NTC samples specified)" > f1 }
  else{ tr=0; for(i=1;i<=ns2;i++){s=sa[i]; tr+=ntc_rem[s]; printf "  %s: %d reads removed\n", s, ntc_rem[s]+0 > f1} printf "  TOTAL removed by NTC subtraction: %d reads\n", tr > f1 }
  close(f1)

  # ===== FILE 2: otuvs_per_sample.txt =====
  f2=P "otuvs_per_sample.txt"
  print "OTUVs per sample" > f2; print "================" > f2; print "" > f2
  for(i=1;i<=ns2;i++){ s=sa[i]; tg=(s in is_ntc)?"  [NTC]":""; dt=(lowdep[s]==1)?"  [LOW DEPTH]":""
    raw=raw_tot[s]+0; ntc_r=ntc_rem[s]+0; low_r=low_rem[s]+0
    # phylum filter loss = post-lowcount total minus phylum-gated OTUV total
    phy_r=(PHYLUM_SET==1)?((pl[s]+0)-(ovg_tot[s]+0)):0; if(phy_r<0)phy_r=0
    # taxonomy(unassigned) loss now only applies within the gated set; with the
    # strict global gate, Unassigned OTUVs are already dropped, so this is 0.
    tax_r=(DO_TAX==1)?((ovg_tot[s]+0)-(otu_tot_s[s]+0)):0; if(tax_r<0)tax_r=0
    if(disp[s]+0==0){ print s tg dt ": All reads from sample were filtered" > f2 }
    else if(DO_TAX==1){
      printf "%s%s%s: %d OTUVs, %d OTUs\n", s, tg, dt, ov_n[s]+0, otu_n[s]+0 > f2
      for(k in otu_r){ split(k,kk,SUBSEP); if(kk[1]==s){ otu=kk[2]
        printf "    - %s (%s): %d reads, %d OTUVs\n", otu, idcov_range(otu), otu_r[k], otu_ovn[k] > f2 } }
    }
    else{
      printf "%s%s%s: %d OTUVs\n", s, tg, dt, ov_n[s]+0 > f2
      for(k in ov_r){ split(k,kk,SUBSEP); if(kk[1]==s){ o=kk[2]
        idcov=(an_pid[o]!=""&&an_pid[o]!="NA")?sprintf("%.0f%% ID, %.0f%% cov", an_pid[o]+0, 100*(an_cov[o]+0)):"no hit"
        lc=(o in an_LCA && an_LCA[o]!="")?an_LCA[o]:"Unavailable"
        printf "    - %s (%s; LCA: %s): %d reads\n", o, idcov, lc, ov_r[k] > f2 } }
    }
    printf "    Reads removed by NTC filtering:    %d (%s)\n", ntc_r, fpct(ntc_r,raw) > f2
    printf "    Reads removed by low-count filter: %d (%s)\n", low_r, fpct(low_r,raw) > f2
    if(PHYLUM_SET==1) printf "    Reads removed by phylum filter:    %d (%s)\n", phy_r, fpct(phy_r,raw) > f2
    if(DO_TAX==1 && tax_r>0) printf "    Reads removed by taxonomy constraint: %d (%s)\n", tax_r, fpct(tax_r,raw) > f2
    print "" > f2
  }
  close(f2)

  # ===== FILE 3: otu_summary.txt =====
  f3=P "otu_summary.txt"
  print "OTU / OTUV summary" > f3; print "==================" > f3; print "" > f3

  if (DO_TAX==1) {
    print "OTUs and Variants (merged by best-hit species)" > f3
    print "----------------------------------------------" > f3
    no=0; for(o in all_otu) ota[++no]=o
    for(i=2;i<=no;i++){k=ota[i];j=i-1;while(j>=1&&ota[j]>k){ota[j+1]=ota[j];j--}ota[j+1]=k}
    for(c=1;c<=no;c++){ otu=ota[c]
      printf "OTU: %s   [%s]   phylum: %s\n", otu, idcov_range(otu), (otu_phy[otu]!=""?otu_phy[otu]:"(none)") > f3
      # member OTUVs (global, deduped)
      nv=split(otu_glob_members[otu], vv, ",")
      # dedup while preserving order
      delete sawm
      for(v=1;v<=nv;v++){ ot=vv[v]; if(ot in sawm) continue; sawm[ot]=1
        idc=(an_pid[ot]!=""&&an_pid[ot]!="NA")?sprintf("%.1f%% ID, %.1f%% cov", an_pid[ot]+0, 100*(an_cov[ot]+0)):"no hit"
        lc=(ot in an_LCA && an_LCA[ot]!="")?an_LCA[ot]:"Unavailable"
        cf=(ot in an_tag && an_tag[ot]!="")?an_tag[ot]:"no-hit"
        printf "    variant %s | %s | LCA: %s | confidence: %s | total reads: %d\n", ot, idc, lc, cf, ovs_tot[ot]+0 > f3
      }
      print "" > f3
    }
    print "" > f3
  }

  print "Best hit + LCA per OTUV" > f3; print "-----------------------" > f3
  no=0; for(o in all_ov) oa[++no]=o
  for(i=2;i<=no;i++){k=oa[i];j=i-1;while(j>=1&&oa[j]>k){oa[j+1]=oa[j];j--}oa[j+1]=k}
  for(c=1;c<=no;c++){ o=oa[c]
    idcov=(an_pid[o]!=""&&an_pid[o]!="NA")?sprintf("%.1f%% ID, %.1f%% cov", an_pid[o]+0, 100*(an_cov[o]+0)):"no hit"
    lc=(o in an_LCA && an_LCA[o]!="")?an_LCA[o]:"Unavailable"
    printf "%s\n", o > f3
    printf "    best hit: %s (%s) | accession: %s\n", (an_S[o]!=""?an_S[o]:"(no species-level name)"), idcov, (an_acc[o]!=""?an_acc[o]:"NA") > f3
    printf "    lineage: P=%s C=%s O=%s F=%s SF=%s G=%s S=%s\n", an_P[o],an_C[o],an_O[o],an_F[o],an_SF[o],an_G[o],an_S[o] > f3
    printf "    LCA (finest): %s\n", lc > f3
    printf "    hits passing coverage: %s | best-bitscore ties: %s | confidence: %s\n", (o in an_nh?an_nh[o]:0), (o in an_nt?an_nt[o]:0), (o in an_tag&&an_tag[o]!=""?an_tag[o]:"no-hit") > f3
    printf "    total reads (all samples): %d\n", ovs_tot[o]+0 > f3
    print "" > f3
  }
  print "Contaminant screening" > f3; print "---------------------" > f3
  if(NTC_LIST==""){ print "(no NTC samples specified)" > f3 }
  else{
    hl="OTUV\ttotal_reads\tntc_reads_total\tntc_fraction"; for(k=1;k<=m;k++)hl=hl OFS ord[k]; print hl > f3
    nc=0; for(o in ovs_tot) if(ovs_tot[o]>0 && (o in all_ov)) ca[++nc]=o
    for(i=2;i<=nc;i++){k=ca[i];kf=ovs_frac[k];j=i-1;while(j>=1&&ovs_frac[ca[j]]<kf){ca[j+1]=ca[j];j--}ca[j+1]=k}
    for(c=1;c<=nc;c++){ o=ca[c]; line=o OFS (ovs_tot[o]+0) OFS (ovs_ntct[o]+0) OFS sprintf("%.3f",ovs_frac[o])
      for(k=1;k<=m;k++) line=line OFS (ovs_ntc[o SUBSEP ord[k]]+0); print line > f3 }
  }
  close(f3)

  # ===== FILE 4: otuv_counts.tsv (wide) =====
  f4=P "otuv_counts.tsv"
  no=0; for(o in all_ov) wa[++no]=o
  for(i=2;i<=no;i++){k=wa[i];j=i-1;while(j>=1&&wa[j]>k){wa[j+1]=wa[j];j--}wa[j+1]=k}
  line="Sample"; for(c=1;c<=no;c++) line=line OFS wa[c]; print line > f4
  for(i=1;i<=ns2;i++){ s=sa[i]; if(s in is_ntc)continue; line=s
    for(c=1;c<=no;c++){ key=s SUBSEP wa[c]; line=line OFS ((key in ov_r && ov_r[key]>0)?ov_r[key]:"NA") } print line > f4 }
  close(f4)

  # ===== FILE 4b: otu_counts.tsv (wide, taxonomy only) =====
  if (DO_TAX==1) {
    f4b=P "otu_counts.tsv"
    no=0; for(o in all_otu) oca[++no]=o
    for(i=2;i<=no;i++){k=oca[i];j=i-1;while(j>=1&&oca[j]>k){oca[j+1]=oca[j];j--}oca[j+1]=k}
    line="Sample"; for(c=1;c<=no;c++) line=line OFS oca[c]; print line > f4b
    for(i=1;i<=ns2;i++){ s=sa[i]; if(s in is_ntc)continue; line=s
      for(c=1;c<=no;c++){ key=s SUBSEP oca[c]; line=line OFS ((key in otu_r && otu_r[key]>0)?otu_r[key]:"NA") } print line > f4b }
    close(f4b)
  }

  # ===== FILE 5: otuv_long.tsv =====
  f5=P "otuv_long.tsv"
  print "Sample","OTUV","Reads","pident","qcov","phylum","class","order","family","subfamily","genus","species","LCA","n_hits","n_ties","confidence" > f5
  for(i=1;i<=ns2;i++){ s=sa[i]
    for(k in ov_r){ split(k,kk,SUBSEP); if(kk[1]==s){ o=kk[2]
      pid=(an_pid[o]!=""&&an_pid[o]!="NA")?sprintf("%.1f",an_pid[o]+0):""
      cov=(an_cov[o]!=""&&an_cov[o]!="NA")?sprintf("%.1f",100*(an_cov[o]+0)):""
      lc=(o in an_LCA && an_LCA[o]!="")?an_LCA[o]:"Unavailable"
      nh=(o in an_nh)?an_nh[o]:0; nt=(o in an_nt)?an_nt[o]:0
      tag=(o in an_tag && an_tag[o]!="")?an_tag[o]:"no-hit"
      print s, o, ov_r[k], pid, cov, an_P[o],an_C[o],an_O[o],an_F[o],an_SF[o],an_G[o],an_S[o], lc, nh, nt, tag > f5 } }
  }
  close(f5)

}
function fpct(part,whole){ if(whole<=0)return "0.0%"; return sprintf("%.1f%%",100.0*part/whole) }
# min-max %ID / %cov range string for a merged OTU
function idcov_range(otu,   ids,covs){
    if (otu in otu_pidmin) {
        if (otu_pidmin[otu]==otu_pidmax[otu]) ids=sprintf("%.1f%% ID", otu_pidmin[otu])
        else ids=sprintf("%.1f-%.1f%% ID", otu_pidmin[otu], otu_pidmax[otu])
    } else ids="no hit"
    if (otu in otu_covmin) {
        if (otu_covmin[otu]==otu_covmax[otu]) covs=sprintf("%.1f%% cov", otu_covmin[otu])
        else covs=sprintf("%.1f-%.1f%% cov", otu_covmin[otu], otu_covmax[otu])
        return ids ", " covs
    }
    return ids
}
'
log "Done. Reports written with prefix: \"${P}\""