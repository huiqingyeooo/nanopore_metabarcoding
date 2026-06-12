################################
### Set paths and variables
################################

PROJECT_NAME=filarioidDiv4
FASTQ_NAME=barcode44

WORKING_DIR=/work/soghigian_lab/huiqing.yeo/metabarcoding
UMI_DEMULTIPLEX_FILE=umi_demultiplex_CLepFol #ensure that file in folder ends with suffix .csv
#SAMPLE_LIST=list_coi.txt
SAMPLE_INFO=${PROJECT_NAME}_sample_info.csv

SEQUENCING_RUN=2026-05-08_metabarcoding
SCRIPTS=/work/soghigian_lab/huiqing.yeo/metabarcoding/scripts

QUALITY=18 #Phred score quality
NEG_CONTROL= #sample name of negative control
COUNT_THRESHOLD= #specify read counts below threshold %
DATABASE_DIR=/work/soghigian_lab/huiqing.yeo/metabarcoding/ref_database/coi
TAXDUMP_DIR=/work/soghigian_lab/huiqing.yeo/metabarcoding/ref_database/
#DB_NAME=nbp-simuliids
#DB_NAME=nem_cerv_hs_dip 
EMU_PERC_ID=97
EMU_MIN_LEN=100

### Setting up paths
METABARCODING_FOLDER=${WORKING_DIR}/${PROJECT_NAME}
FASTQ_DIR=${WORKING_DIR}/raw_fastq_files/${SEQUENCING_RUN}
DEMULTIPLEX_DIR=${METABARCODING_FOLDER}/0_fastq_files
TRIMMED_FASTQ_DIR=${METABARCODING_FOLDER}/1_fastq_trimmed
QC_DIR=${METABARCODING_FOLDER}/2_QC
EMU_DIR=${METABARCODING_FOLDER}/3_emu
VSEARCH_DIR=${METABARCODING_FOLDER}/4_vsearch
VSEARCH_CLUSTERING_THRESHOLD=97
GENBANK_DATABASE=/work/soghigian_lab/databases/gb_nt/nt_euk

#########################################
### concatenate fastq files
#########################################
salloc --mem=150G -c 1 -N 1 -n 8 -t 05:00:00

cd ${WORKING_DIR}
mkdir raw_fastq_files/${SEQUENCING_RUN}
cd raw_fastq_files/${SEQUENCING_RUN}

for i in {42..47};
do
echo barcode$i
cat /work/soghigian_lab/data/nanopore/20260506_metabarcoding_2026-05-08_20-28-46/barcode$i.fastq.gz \
/work/soghigian_lab/data/nanopore/20260508_metabarcoding_2_2026-05-08_20-29-51/barcode$i.fastq.gz > \
barcode$i.fastq.gz
done

# set up folders
mkdir ${METABARCODING_FOLDER}
cd ${METABARCODING_FOLDER}
mkdir 0_fastq_files 1_fastq_trimmed 2_QC 3_emu 4_vsearch

# upload umi demultiplexing_file.csv and sample_info.csv into 0_fastq_files directory
# upload sample list into 0_fastq_files directory too
# note: sample list can be a subset of samples listed in sample_info.
# for e.g., a run has samples with both COI and ITS2 genes, there will be two sample lists
# another e.g., run has samples which requires matching to two separate databases

#########################################
### umi demultiplexing
#########################################
cd ${SCRIPTS}
# usage: sbatch 1_umi_demultiplex.slurm --output path/to/working/dir/slurm-%j.out \
# <fastq_prefix> <umi_file_name> <fastq_dir> <working_dir>

sbatch --output ${DEMULTIPLEX_DIR}/slurm-%j.out 01_umi_demultiplex.slurm \
${FASTQ_NAME} ${UMI_DEMULTIPLEX_FILE} ${FASTQ_DIR} ${DEMULTIPLEX_DIR}

### create sample lists for demultiplexing
# note: each list is specific to each gene

#########################################
### trim primers
#########################################
env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
sbatch --output ${TRIMMED_FASTQ_DIR}/slurm-%j.out \
02_filter_trimPrimers.slurm :::: ${DEMULTIPLEX_DIR}/list_coi.txt ::: ${QUALITY} ::: COI_INSECT ::: ${DEMULTIPLEX_DIR} ::: ${TRIMMED_FASTQ_DIR}

# check log files
for i in $(ls ${TRIMMED_FASTQ_DIR}/slurm-*.out);
do 
col1=`echo ${i}`
col2=`tail -n 1 ${i}`
echo $col1 $col2 >> ${TRIMMED_FASTQ_DIR}/slurm_summary.txt
done

#########################################
### QC check
#########################################
env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
sbatch --output ${QC_DIR}/slurm_output/slurm-%j.out \
03_nanoplot_QC.slurm :::: ${DEMULTIPLEX_DIR}/list_coi.txt ::: ${QUALITY} ::: ${DEMULTIPLEX_DIR} ::: ${TRIMMED_FASTQ_DIR} ::: ${QC_DIR}

#########################################
### summarize nanoplot QC results
#########################################
sbatch --output ${QC_DIR}/slurm_output/slurm-%j.out \
04_summarise_QC.slurm ${QC_DIR}

##############################################################################################
### Identify samples with read count below negative control and below a specified threshold
##############################################################################################
# this generates a list of samples with reads above negative control (sample_read_count_above_nc.txt)
# and generates list of samples with reads above specified threshold (sample_read_count_above_cutoff.txt)
# list can be used as input for the next step

# Did not filter samples in this run
sbatch --output ${WORKING_DIR}/slurm-%j.out \
05_filter_samples.slurm ${NEG_CONTROL} q${QUALITY} ${COUNT_THRESHOLD} ${TRIMMED_FASTQ_DIR}

#########################################
### run emu
#########################################
# batch code

# To generate sample list from files in a folder
SAMPLE_LIST=$(ls ${TRIMMED_FASTQ_DIR} | grep q"${QUALITY}" | grep "trimmed" | sed "s#_trimmed.fastq.gz##g")

# Run emu on filtered samples
SAMPLE_LIST=$(awk '{print $1}' ${TRIMMED_FASTQ_DIR}/sample_read_count_above_cutoff.txt | grep q"${QUALITY}" | grep "trimmed" | sed "s#_trimmed.fastq.gz##g")

cd ${SCRIPTS}
# usage: parallel sbatch --output ${WORKING_DIR}/slurm-%j.out 6_emu.slurm <db_name> <percentage_id_cutoff> <min_length_cutoff> <sample> \
# </path/to/database> </path/to/fastq/files> </path/for/emu/output/files>

env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
sbatch --output ${EMU_DIR}/slurm_output/slurm-%j.out \
06_emu.slurm ::: ${DB_NAME} ::: ${EMU_PERC_ID} ::: ${EMU_MIN_LEN} :::: ${SAMPLE_LIST} ::: ${DATABASE_DIR} ::: ${TRIMMED_FASTQ_DIR} ::: ${EMU_DIR}


# run emu using sample list in a text file
SAMPLE_LIST=list_coi_CBAY0063.txt
DB_NAME=nbp-simuliids

SAMPLE_LIST=list_coi_other.txt
DB_NAME=nem_cerv_hs_dip

sed "s/$/_q$QUALITY/" ${DEMULTIPLEX_DIR}/${SAMPLE_LIST} > ${DEMULTIPLEX_DIR}/${SAMPLE_LIST}2

env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
sbatch --output ${EMU_DIR}/slurm_output/slurm-%j.out \
06_emu.slurm ::: ${DB_NAME} ::: ${EMU_PERC_ID} ::: ${EMU_MIN_LEN} :::: ${DEMULTIPLEX_DIR}/${SAMPLE_LIST}2 ::: ${DATABASE_DIR} ::: ${TRIMMED_FASTQ_DIR} ::: ${EMU_DIR}


# check log files
for i in $(ls ${EMU_DIR}/slurm_output/slurm-*.out);
do 
col1=`echo ${i}`
col2=`tail -n 1 ${i}`
echo $col1 $col2 >> ${EMU_DIR}/slurm_output/slurm_summary.txt
done

###########################################
### concatenate and emu summarise results
###########################################
#usage: sbatch --output ${EMU_DIR}/slurm-%j.out \
#7_concat_tsv.slurm <percentage_id_cutoff> <min_length_cutoff> </path/to/emu/output/file> 

sbatch --output ${EMU_DIR}/slurm_output/slurm-%j.out \
07_concat_tsv.slurm ${EMU_PERC_ID} ${EMU_MIN_LEN} ${EMU_DIR} ${DEMULTIPLEX_DIR} ${SAMPLE_INFO}


#############################################
### investigate filtered and unmapped reads
#############################################
mkdir ${VSEARCH_DIR}

### first rename headers in fastq files with sample names
# batch code
# Run on all trimmed samples
SAMPLE_LIST=$(ls ${TRIMMED_FASTQ_DIR} | grep "${QUALITY}" | grep "trimmed" | sed "s#_q18_trimmed.fastq.gz##g")

# Run on filtered samples
SAMPLE_LIST=$(awk '{print $1}' ${TRIMMED_FASTQ_DIR}/sample_read_count_above_cutoff.txt | grep "${QUALITY}" | grep "trimmed" | sed "s#_q18_trimmed.fastq.gz##g")

env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
sbatch --output ${VSEARCH_DIR}/slurm_output/slurm-%j.out \
08_relabel.slurm ::: ${SAMPLE_LIST} ::: q${QUALITY} ::: ${EMU_DIR}/pid97_len100 ::: ${VSEARCH_DIR}

### Cluster sequences
sbatch --output ${VSEARCH_DIR}/slurm-%j.out \
09_cluster.slurm ${VSEARCH_DIR} ${PROJECT_NAME} ${VSEARCH_CLUSTERING_THRESHOLD}

### Blast OTUs
sbatch --output ${VSEARCH_DIR}/slurm-%j.out \
10_blast.slurm ${PROJECT_NAME} ${VSEARCH_CLUSTERING_THRESHOLD} ${VSEARCH_DIR} ${GENBANK_DATABASE} ${TAXDUMP_DIR}
