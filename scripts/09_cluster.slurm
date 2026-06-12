#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=50G
#SBATCH --output=slurm-cluster-%j.out

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/vsearch

VSEARCH_DIR=${1}
PROJECT_NAME=${2}
CLUSTERING_THRESHOLD=${3}

####### script ##################################
cd ${VSEARCH_DIR}

cat *_relabelled.fastq > ${PROJECT_NAME}.fastq
#rm *_relabelled*

# dereplicate 
vsearch --fastx_uniques ${PROJECT_NAME}.fastq \
--sizein \
--sizeout \
--fasta_width 0 \
--uc ${PROJECT_NAME}_derep.uc \
--fastqout ${PROJECT_NAME}_derep.fastq

#--minuniquesize 2 \

echo "vsearch dereplicate done"

# preclustering for chimera detection
vsearch --cluster_size ${PROJECT_NAME}_derep.fastq \
--threads 8 \
--id 0.98 \
--strand plus \
--sizein \
--sizeout \
--clusterout_sort \
--fasta_width 0 \
--uc ${PROJECT_NAME}_derep_preclustered.uc \
--clusterout_id \
--centroids ${PROJECT_NAME}_derep_preclustered.fasta

echo "preclustering for chimera detection done"

# de novo chimera detection
vsearch --uchime_denovo ${PROJECT_NAME}_derep_preclustered.fasta \
--sizein \
--sizeout \
--fasta_width 0 \
--nonchimeras ${PROJECT_NAME}_derep_preclustered_nonchimeras.fasta

echo "de novo chimera detection done"

# Cluster at specified threshold, and relabel with OTU_n, generate OTU table

vsearch --cluster_size ${PROJECT_NAME}_derep_preclustered_nonchimeras.fasta \
--threads 8 \
--id 0.${CLUSTERING_THRESHOLD} \
--strand plus \
--sizein \
--sizeout \
--fasta_width 0 \
--relabel OTU_ \
--uc ${PROJECT_NAME}_derep_${CLUSTERING_THRESHOLD}otu.uc \
--centroids ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otu.fasta \
--otutabout ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab.txt

echo "clustering at ${CLUSTERING_THRESHOLD}% done"

### summarize OTU table

# summarize counts across columns
awk -F'\t' '
BEGIN { OFS="\t" }

NR == 1 {
    # Build column → sample mapping
    for (i = 2; i <= NF; i++) {
        split($i, a, "_")
        sample = a[1]
        if (!(sample in idx)) {
            idx[sample] = ++n
            order[n] = sample
        }
        col[i] = idx[sample]
    }

    # Print new header
    header = "OTU"
    for (i = 1; i <= n; i++) header = header OFS order[i]
    print header
    next
}

{
    # Initialize sums
    for (i = 1; i <= n; i++) sum[i] = 0

    # Sum replicate columns
    for (i = 2; i <= NF; i++)
        sum[col[i]] += $i

    # Print summarized row
    line = $1
    for (i = 1; i <= n; i++) line = line OFS sum[i]
    print line
}
' ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab.txt > ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized.txt

echo "Summarizing of OTU table done"

# reorder OTUs. e.g. OTU_1, OTU_2, OTU_3
rm ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized_sorted.txt
head -n 1 ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized.txt > ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized_sorted.txt
tail -n +2 ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized.txt | sort -t $'\t' -k1,1V >> \
${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized_sorted.txt

echo "Reordering OTU table done"

rm ${PROJECT_NAME}_derep.*
rm ${PROJECT_NAME}_${CLUSTERING_THRESHOLD}otutab_summarized.txt

echo "Job finised with exit code $? at: `date`"
