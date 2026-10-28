#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=120G

####### Notes ######################
# 21 October 2026
# Base script from John Soghigian
# I am running this on interactively on ARC to test it out. I am using the files for the 20260923_COI_MB017-MB020 nanopore run.

####### Description of Pipeline #######
# At present, this pipeline is essentially a plate_1 loop script. It can be run in different ways, but the main goal is for a job per plate to be run. 

# Note that the pipeline assumes you have already done a few things. Notably:
# 1) Used fastplong or another trimming software to output reads of decent quality.  We use Q15 or Q18.
# 2) Prepared input TSVs that contain information on the sample, namely sample ID and UMIs.
# In both cases, a companion script has been provided to help.

#Example of the command to generate the trimmed reads
#source ~/software/init-conda-sl
#conda activate /work/soghigian_lab/apps/conda/envs/chopper_v0.11.0

#FASTQ_DIR=/work/soghigian_lab/data/20260923_NBD_library1_library2_fastq_2026-09-25_14-21-19/20260923_library1_HY_OC/
#for i in 50 51 52 53 55;
#do
#echo barcode${i}
#pigz -dc ${FASTQ_DIR}/barcode${i}/barcode${i}.fastq.gz | chopper -q 15 --minlength 700 --maxlength 900 | pigz > barcode${i}_q15.fastq.gz
#done

# Submit jobs with parallel 
# 
#parallel --link --delay 1 sbatch pipeline_polish_cluster.sh ::: MB017 MB018 MB019 MB020 NTC ::: barcode50 barcode51 barcode52 barcode53 barcode55

####### Set environment variables ###############

### Change these for your particular project and job:
FASTQ_FILE=/work/soghigian_lab/huiqing.yeo/metabarcoding/20260923_COI_MB017-MB020/pipeline_v2/${2}_q15.fastq.gz
plateid=${1}
proj=plate_${plateid}
THREADS=8
phylum=arthropoda
demuxfile=/work/soghigian_lab/huiqing.yeo/metabarcoding/20260923_COI_MB017-MB020/pipeline_v2/plates_umi_demux_${plateid}.csv

SCRIPTSDIR=/work/soghigian_lab/huiqing.yeo/metabarcoding/scripts/nanopore_metabarcoding/pipeline02_polish_cluster
F_PRIMER=RKTCAACMAATCATAAAGATATTGG #Forward primer, 5' -> 3' direction
R_PRIMER=TAAACTTCWGGRTGWCCAAAAAATCA #Reverse primer, 5' -> 3' direction
WORKING_DIR=/work/soghigian_lab/huiqing.yeo/metabarcoding/20260923_COI_MB017-MB020/pipeline_v2

cd ${WORKING_DIR}
### Define the NTCs - Comment out this line if you'd like to NOT use NTCs.

#ntcs=$(grep "Neg" plate_1.tsv | awk -F"\t" '{ print $1 }' | tr '\n' ',')
#ntcs=plateTEMPLATE_95_PTEMPLATE_Neg1,plateTEMPLATE_96_PTEMPLATE_Neg2

### These are internal settings that will change some output
#A couple of quick terms to know - in this pipeline, because we are using noisy Nanopore reads, we do not generate ASVs specifically but what we are calling "OTUVs", which are fine-grained read clusters; they are similar to ASVs, but they are generated from read clustering in isONclust, then polished.
#This is the number of reads used at most for generating consensus sequences from a contig
MAXDRAFT=500
#This defines the minimum number of reads per cluster. If an analysis produces too many OTUVs, consider increasing this.
minclust=50

#This is a de-noising optional filter that will adjust final reported read counts per OTUV.
lowcount_num=50

#This just provides a warning for samples below a certain read depth
warning=500

#Required databases for the script
acc2taxid=/work/soghigian_lab/databases/gb_nt/acc2taxid.new
taxonkit=/work/soghigian_lab/apps/taxonkit
taxdump=/work/soghigian_lab/databases/taxdump

########## PIPELINE STARTS BELOW ##########
#Change to YOUR conda environment!
source ~/software/init-conda-sl

### Load required programs
module load gcc
conda activate /work/soghigian_lab/apps/conda/envs/vsearch

SEQTK=/work/soghigian_lab/apps/seqtk/seqtk
VSEARCH=$(which vsearch)
isONclust=/work/soghigian_lab/apps/isONclust3/target/release/isONclust3
SPOA=/work/soghigian_lab/apps/spoa/build/bin/spoa
SEQKIT=/work/soghigian_lab/apps/seqkit
MINIMAP=/work/soghigian_lab/apps/minimap2-2.31_x64-linux/minimap2
FASTPLONG=/work/soghigian_lab/apps/fastplong/fastplong
conda activate /work/soghigian_lab/apps/conda/envs/cutadaptenv

### Raise the number of files:
ulimit -S -n 4096

mkdir -p ${WORKING_DIR}/${proj}
cd ${WORKING_DIR}/${proj}

# make_plate_fastas.py script: takes UMIs and primers and generates a fasta file with linked UMI+F and reverse complimented R+UMI
# the first 8 nucleotides of both forward and reverse primers are used as anchors to ensure that the UMIs are not matched inside the amplicon.
python3 ${SCRIPTSDIR}/make_plate_fastas.py $demuxfile -o fasta --F ${F_PRIMER} --R ${R_PRIMER}

bash ${SCRIPTSDIR}/demux_plate.sh $plateid $FASTQ_FILE fasta demux 12 1500 1 $F_PRIMER $R_PRIMER  # repeat per plate

conda activate /work/soghigian_lab/apps/conda/envs/vsearch
echo merging reads...

### Cluster reads - begin by merging sample reads:
cat demux/final/*.fastq.gz > all.fq.gz
gunzip all.fq.gz
#rm $readpath/all.fastq.gz
echo Clustering reads...

$isONclust --fastq all.fq --mode ont -n $minclust --min-shared-minis 0.9 --outfolder ${WORKING_DIR}/${proj}

echo Generating consensus from clusters...

cd ${WORKING_DIR}/${proj}

mkdir -p polished

for f in clustering/fastq_files/*.fastq; do
    id=$(basename "$f" .fastq)

    # count reads in this cluster
    n=$($SEQKIT stats -T "$f" | awk 'NR==2{print $4}')

    if [ "$n" -gt "$MAXDRAFT" ]; then
        echo "[$id] $n reads (>$MAXDRAFT): subsetting to $MAXDRAFT for draft"
        $SEQKIT sample -n "$MAXDRAFT" "$f" > "polished/${id}_sub.fastq"
    else
        echo "[$id] $n reads (<=$MAXDRAFT): using all reads for draft"
        cp "$f" "polished/${id}_sub.fastq"
    fi

    # draft consensus from the (possibly subsetted) reads
    $SPOA "polished/${id}_sub.fastq" > "polished/${id}_draft.fasta"

    # relabel header and collect
     sed "s/^>.*/>cluster_${id}/" "polished/${id}_draft.fasta" >> consensus.fasta
    rm -f "polished/${id}_sub.fastq"   # tidy up the temp subset
done

echo Polishing cluster contigs

$MINIMAP -t $THREADS -x map-ont --secondary=no \
  consensus.fasta \
  clustering/sorted.fastq > polished/all_to_reps.paf
 
conda activate /work/soghigian_lab/apps/conda/envs/medaka_2.2.2

medaka_consensus -t $THREADS \
  -i clustering/sorted.fastq \
  -d consensus.fasta \
  -o medaka_consensusnoR_out \
  -m r1041_e82_400bps_sup_v5.2.0

conda activate /work/soghigian_lab/apps/conda/envs/vsearch

#Because of polishing, we could end up with some identical sequences. Importantly, dereplication won't work on these, because they are often not totally overlapping. So instead, we cluster... with a 100% identity threshold.

log Clustering polished contigs at 100% identity to remove duplicates

mkdir prelim_otus

vsearch --cluster_size consensus.fasta \
  --id 1 \
  --centroids prelim_otus/otus_unique.fasta \
  --uc prelim_otus/clusters.unique \
  --sizein --sizeout

# now let's remove chimeras. Save the chimeras as
vsearch --uchime3_denovo prelim_otus/otus_unique.fasta \
  --nonchimeras otuvs.fasta \
  --chimeras prelim_otus/otus_chimeras.fasta

###Now let's examine read depth across OTUVs...

log "Counting reads per OTUV per sample"

mkdir -p counts

for FQ in ./demux/final/*.fastq.gz; do
  RUN=$(basename "$FQ" .fastq.gz)
  $MINIMAP -t $THREADS -x map-ont --secondary=no \
    otuvs.fasta "$FQ" \
    > counts/${RUN}.paf
done

{
  echo -e "Sample\tOTUV\tReads"
  for P in counts/*.paf; do
    RUN=$(basename "$P" .paf)
    awk 'BEGIN{FS="\t"; OFS="\t"}
         {key=$6"\t"$1; if(!(key in seen)) { seen[key]=1; c[$6]++ }}
         END{for(h in c) print "'"$RUN"'", h, c[h]}' "$P"
  done
} > otuv_counts.tsv

#rm -r medaka_consensusnoR_out/
#rm -r clustering/
#rm -r polished/

### Now we will determine the identity of our sequneces.  For our purposes, we'll use blast. There are faster alternatives, but, this is comprehensive, and provides a taxonomy ID we can retrieve.
log "Running BLAST to identify OTUVs"

BLASTDB=/work/soghigian_lab/databases/gb_nt
PATH=$PATH:/work/soghigian_lab/apps/ncbi-blast-2.13.0+/bin
blastn -query otuvs.fasta -num_threads 16 \
  -db /work/soghigian_lab/databases/gb_nt/nt_euk_prok_new \
  -task megablast \
  -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
  -out haplotypes_vs_COIrefs.b6

log "Summarizing run results"
sh $SCRIPTSDIR/summarize_run_2.1.sh -c otuv_counts.tsv \
  -b haplotypes_vs_COIrefs.b6 --acc2taxid=$acc2taxid \
  --taxonkit=$taxonkit \
  --taxdump=$taxdump \
  --phylum=$phylum \
  --ntc=$ntcs \
  --lowcount=$lowcount_num --mindepth=$warning \
  --taxonomy -o ${proj}_

#rm ../all.fq