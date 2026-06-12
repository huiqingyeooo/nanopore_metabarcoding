#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/chopper_v0.11.0

SAMPLE=$1
QUALITY=$2
MARKER=$3
WORKING_DIR=$4
TRIMMED_FASTQ_DIR=$5

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 \
#sbatch --output ${WORKING_DIR}/slurm-%j.out \
#2_filter_trimPrimers.slurm :::: ${WORKING_DIR}/list_coi.txt ::: 18 ::: COI_NEMATODE ::: ${WORKING_DIR}

####### script ##################################
# reference script from:
# https://github.com/henrimdias/emu-microbiome-HPC/blob/main/emu_pipeline/scripts/run_nanotools_pipeline.slurm
cd ${WORKING_DIR}

if [[ "$MARKER" == "COI_NEMATODE" ]]; then
    # filter reads
    chopper -q ${QUALITY} --minlength 540 --maxlength 740 -i ${SAMPLE}.fastq.gz | gzip > ${SAMPLE}_q${QUALITY}.fastq.gz
    # trim primers
    conda activate /work/soghigian_lab/apps/conda/envs/cutadaptenv
    cutadapt --cores=4 -e 4 -g TGATTGGTGGTTTTGGTAA...GATATTGATACTCGTACTTAT \
    -o ${TRIMMED_FASTQ_DIR}/${SAMPLE}_q${QUALITY}_trimmed.fastq.gz --revcomp --discard-untrimmed --overlap 18 \
    ${SAMPLE}_q${QUALITY}.fastq.gz

elif [[ "$MARKER" == "COI_INSECT" ]]; then
    # filter reads
    chopper -q ${QUALITY} --minlength 540 --maxlength 740 -i ${SAMPLE}.fastq.gz | gzip > ${SAMPLE}_q${QUALITY}.fastq.gz
    # trim primers
    conda activate /work/soghigian_lab/apps/conda/envs/cutadaptenv
    cutadapt --cores=4 -e 4 -g RKTCAACMAATCATAAAGATATTGG...TGRTTTYTTGGTCACCCTGAAGTTTA \
    -o ${TRIMMED_FASTQ_DIR}/${SAMPLE}_q${QUALITY}_trimmed.fastq.gz --revcomp --discard-untrimmed --overlap 21 \
    ${SAMPLE}_q${QUALITY}.fastq.gz
    
elif [[ "$MARKER" == "ITS2_NEMATODE" ]]; then
    # filter reads
    chopper -q ${QUALITY} --minlength 900 --maxlength 1500 -i ${SAMPLE}.fastq.gz | gzip > ${SAMPLE}_q${QUALITY}.fastq.gz
    # trim primers
    conda activate /work/soghigian_lab/apps/conda/envs/cutadaptenv
    cutadapt --cores=4 -e 4 -g GTAGGTGAACCTGCGGAAGGATCATT...AGCGGAGGAAAAGAAACTAA \
    -o ${TRIMMED_FASTQ_DIR}/${SAMPLE}_q${QUALITY}_trimmed.fastq.gz --revcomp --discard-untrimmed --overlap 21 \
    ${SAMPLE}_q${QUALITY}.fastq.gz

else
    echo "Error: marker must be COI_INSECT, COI_NEMATODE or ITS2_NEMATODE, got '${MARKER}'"
    exit 1
fi
echo "Job finised with exit code $? at: `date`"
