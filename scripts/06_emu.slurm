#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem=20G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/emu-3.6.2
EMU=/work/soghigian_lab/apps/conda/envs/emu-3.6.2/emu

DB_NAME=${1}
PID=${2}
LEN=${3}
SAMPLE=${4}
DATABASE=${5}
FASTQ_FILES=${6}
OUTPUT_DIR=${7}/pid${2}_len${3}

echo "starting run at: `date`"
####### script ##################################
mkdir ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# fasta file
#emu abundance ${sample}.fasta \
#--db ${database}/${db_name} --threads 24 --output-dir ${sample} \
#--type map-ont --keep-counts --keep-read-assignments

# fastq file
${EMU} abundance --type sr ${FASTQ_FILES}/${SAMPLE}_trimmed.fastq.gz \
--db ${DATABASE}/${DB_NAME} --threads 12 --output-dir ${OUTPUT_DIR}/${SAMPLE} \
--min-pid ${PID} --min-align-len ${LEN} \
--type map-ont --keep-counts --keep-read-assignments \
--output-unclassified

echo "Job finised with exit code $? at: `date`"
