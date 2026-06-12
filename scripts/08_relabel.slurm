#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --output=slurm-%j.out

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/vsearch

SAMPLE=${1}
QUALITY=${2}
EMU_DIR=${3}
VSEARC_DIR=${4}

####### script ##################################
cd ${VSEARC_DIR}
echo "relabelling fastq sequences for ${SAMPLE}"

# relabel headers of fastq files to include sample names
# ensures preservation of sample info to each sequence when concatenating fastq files
# fastq_max 93 as any number above that throws an error

vsearch --fastx_filter ${EMU_DIR}/${SAMPLE}_${QUALITY}/${SAMPLE}_${QUALITY}_trimmed.fastq_filtered_mapped.fastq \
--fastqout ${SAMPLE}_${QUALITY}_trimmed.fastq_filtered_mapped_relabelled.fastq \
--relabel_keep --relabel ${SAMPLE}filteredmapped_ --fastq_qmax 93

vsearch --fastx_filter ${EMU_DIR}/${SAMPLE}_${QUALITY}/${SAMPLE}_${QUALITY}_trimmed.fastq_unclassified_mapped.fastq \
--fastqout ${SAMPLE}_${QUALITY}_trimmed.fastq_unclassified_mapped_relabelled.fastq \
--relabel_keep --relabel ${SAMPLE}unclassifiedmapped_ --fastq_qmax 93

vsearch --fastx_filter ${EMU_DIR}/${SAMPLE}_${QUALITY}/${SAMPLE}_${QUALITY}_trimmed.fastq_unmapped.fastq \
--fastqout ${SAMPLE}_${QUALITY}_trimmed.fastq_unmapped_relabelled.fastq \
--relabel_keep --relabel ${SAMPLE}unmapped_ --fastq_qmax 93

echo "Job finised with exit code $? at: `date`"
