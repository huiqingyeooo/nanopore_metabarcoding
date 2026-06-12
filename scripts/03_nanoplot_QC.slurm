#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/nanoplot

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
#sbatch --output ${QC_DIR}/slurm-%j.out \
#03_nanoplot_QC.slurm :::: ${DEMULTIPLEX_DIR}/list_coi.txt ::: ${QUALITY} ::: ${DEMULTIPLEX_DIR} ::: ${TRIMMED_FASTQ_DIR} ::: ${QC_DIR}

SAMPLE=${1}
QUALITY=${2}
DEMULTIPLEX_DIR=${3}
TRIMMED_FASTQ_DIR=${4}
QC_DIR=${5}

####### script ##################################
# reference script from:
# https://github.com/henrimdias/emu-microbiome-HPC/blob/main/emu_pipeline/scripts/run_nanotools_pipeline.slurm

# QC reads before filtering
NanoPlot --fastq ${DEMULTIPLEX_DIR}/${SAMPLE}.fastq.gz --threads 8 \
--outdir ${QC_DIR}/raw/${SAMPLE} \
--legacy hex dot \
--maxlength 2000 -f pdf

# QC reads after filtering
NanoPlot --fastq ${TRIMMED_FASTQ_DIR}/${SAMPLE}_q${QUALITY}_trimmed.fastq.gz --threads 8 \
--outdir ${QC_DIR}/filtered/${SAMPLE} \
--legacy hex dot \
--maxlength 2000 -f pdf

echo "Job finised with exit code $? at: `date`"
