#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --output=slurm-blast-%j.out

####### Set environment variables ###############
source ~/software/init-conda
#conda activate /work/soghigian_lab/apps/conda/envs/blast
PATH=/home/huiqing.yeo/software/ncbi-blast-2.14.1+/bin:$PATH

PROJECT_NAME=${1}
VSEARCH_CLUSTERING_THRESHOLD=${2}
VSEARCH_DIR=${3}
BLAST_DATABASE=${4}
TAXDUMP_DIR=${5}

echo "starting run at: `date`"

####### Script ###################################
#this step is for creating a custom database
#cd /work/soghigian_lab/huiqing.yeo/filiaroidDiv/test/db
#makeblastdb -in db.fasta -dbtype nucl -out db
cd ${VSEARCH_DIR}

# download nucl_gb.accession2taxid.gz file from NCBI
# use accession column instead of accession.version column in case of unmatched versions
# only need to run this script once, and the lines can be hashed out for subsequent runs

# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# zcat nucl_gb.accession2taxid.gz | cut -f 1,3 | gzip -c > acc2taxid.tsv.gz

# top 5 blast matches 
blastn -db ${BLAST_DATABASE} \
-query ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu.fasta  -out ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast5.out \
-outfmt "6 qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-max_target_seqs 5

# Obtain top blast result every 5 lines
# prints first line every five lines
awk 'NR % 5 == 1' ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast5.out > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1.out

# use taxonkit to match species ID with accession numbers
conda activate /work/soghigian_lab/apps/conda/envs/taxonkit

# extract accession and deduplicate and remove versions to reduce memory usage

# for customized database
# 12 columns
#OTU_1;size=377	OM401548.1	98.949	571	5	1	69	638	1	571	0.0	1020
#OTU_2;size=2	MG677940.1	87.809	648	77	2	1	646	10	657	0.0	758
#OTU_3;size=2	OM401548.1	95.147	577	15	10	69	638	1	571	0.0	898
#OTU_4;size=28	AJ271617.1	93.569	622	38	2	8	627	7	628	0.0	926

#cut -f 2 sample01_blast_${VSEARCH_CLUSTERING_THRESHOLD}otus.out | csvtk uniq -Ht | csvtk replace -Ht -p '\.\d+$' > acc.txt

# for nt_euk database
# 13 columns (taxid is in 3rd column)
#OTU_1;size=372	gi|1335850942|gb|MG677940.1|	153009	87.809	648	77	2	1	646	10	657	0.0	758
#OTU_1;size=372	gi|2020302796|gb|MW835250.1|	153009	87.809	648	77	2	1	646	13	660	0.0	758
#OTU_1;size=372	gi|149384889|gb|EF521410.1|	153009	87.809	648	77	2	1	646	1	648	0.0	758
#OTU_1;size=372	gi|1032530272|gb|KX132091.1|	153009	87.893	636	75	2	13	646	1	636	0.0	747

cut -f 2 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1.out | awk -F'|' '{print $4}' | csvtk uniq -Ht | csvtk replace -Ht -p '\.\d+$' > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc.out

# grep taxids from acc2taxid.tsv.gz with accession list from blast results (acc.out)
zcat ${TAXDUMP_DIR}/acc2taxid.tsv.gz | grep -w -f ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc.out >  ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc2taxid.out

# download and uncompress ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz, and copy "names.dmp", "nodes.dmp", "delnodes.dmp", and "merged.dmp"
# retrieve species name for taxids

cut -f 2 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc2taxid.out | taxonkit reformat -f '{s}' -I 1 --data-dir ${TAXDUMP_DIR}/taxdump > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_taxid2name.out

# replace spaces with "_" in species names, for easier downstream sorting
sed -i "s# #_#g" ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_taxid2name.out

# append taxids to accessions, append species names to taxids

# for custom blast database
# header for 12 columns
#csvtk add-header -t --names "counts,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" \
#${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1.out \
#    | csvtk mutate -t -f sseqid -n taxid \
#    | csvtk replace -t -k ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc2taxid.out -f taxid -p '(.+)\.\d+' -r '{kv}' \
#    | csvtk mutate -t -f taxid -n species \
#    | csvtk replace -t -k ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_taxid2name.out -f species -p '(.+)' -r '{kv}' > \
#    ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id.out

#genbank nt_euk database
csvtk add-header -t --names "counts,sseqid,taxid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" \
${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1.out \
    | csvtk replace -t -k ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc2taxid.out -f taxid -p '(.+)\.\d+' -r '{kv}' \
    | csvtk mutate -t -f taxid -n species \
    | csvtk replace -t -k ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_taxid2name.out -f species -p '(.+)' -r '{kv}' > \
    ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id.out
#| csvtk mutate -t -f sseqid -n taxid \

# split first column (OTU_x; size=x)
awk '{sub (/;/, OFS)} 1' OFS="\t" ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id.out > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format.out

### join otutable and blast results (n=1)
tail -n +2 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format.out > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format2.out

tail -n +2 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted.txt > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted_format2.txt

head -n 1 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format.out > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format_header.out

head -n 1 ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted.txt > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted_header.txt
sed -i "s#\t# #g" ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted_header.txt

paste ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted_header.txt ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format_header.out > tmp_header.txt

join ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otutab_summarized_sorted_format2.txt ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_id_format2.out > tmp_wo_header.txt

cat tmp_header.txt tmp_wo_header.txt > ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}_otutable.final.csv

sed -i "s/size=//g" ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}_otutable.final.csv

# housekeeping
rm ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc.out
rm ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_acc2taxid.out
rm ${PROJECT_NAME}_${VSEARCH_CLUSTERING_THRESHOLD}otu_blast1_taxid2name.out
rm *format*
rm *header*
rm *tmp*
echo "Job finised with exit code $? at: `date`"
