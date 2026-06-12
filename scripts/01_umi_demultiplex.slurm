#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=8:00:00
#SBATCH --mem=80G

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: sbatch 1_umi_demultiplex.slurm <fastq_prefix> <umi_file_name> <fastq_dir> <working_dir>"
    echo "  e.g. sbatch 1_umi_demultiplex.slurm barcode01 umi_demultiplex_barcode01 /path/to/fastq/files/dir /path/to/working/dir"
    exit 1
fi

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/cutadaptenv

BARCODE=$1
UMI_FILE=$2
FASTQ_DIR=$3
WORKING_DIR=$4

FASTQ_FOLDER=${FASTQ_DIR}/${BARCODE}.fastq.gz

####### script ##################################
# single sample demultiplex
#name=sample5
#cutadapt --cores=4 -e 2 -g ACACATCTCGTGAGAG...GCAGACTCTCACACGC \
#--discard-untrimmed -o ${name}.fastq.gz --revcomp --overlap 16 \
#barcode02.fastq.gz

# prepare umi_demultipex.fasta
python convert_umi_fasta.py ${WORKING_DIR}/${UMI_FILE}.csv ${WORKING_DIR}/${UMI_FILE}.fasta

# batch sample demultiplex
cd ${WORKING_DIR}

cutadapt --cores=24 -e 2 -g "file:${UMI_FILE}.fasta" \
--discard-untrimmed --revcomp --overlap 16 -o {name}.fastq.gz \
${FASTQ_FOLDER}

echo "Job finised with exit code $? at: `date`"
