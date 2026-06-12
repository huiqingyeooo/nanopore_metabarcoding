#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=80G

####### Set environment variables ###############

# Define parameters
NEG_CONTROL=${1}
QUALITY=${2}
CUTOFF=${3}
# 0.0005 = 0.05% of total reads
WORKING_DIR=${4}

####### script ##################################
cd ${WORKING_DIR}

# This script filters out samples with read counts lower than negative control
rm read_counts.txt

# loop through samples to count reads in *_trimmed.fastq.gz file
for i in *_${QUALITY}_trimmed.fastq.gz;
do
col1=`echo ${i}`
col2=`echo $(zcat ${i}|wc -l)/4 | bc`
echo $col1 $col2 >> read_counts.txt
done

# Calculate total number of reads
TOTAL_READS=$(awk '{s+=$2} END {print s}' read_counts.txt)
echo "Total number of reads across all samples: $TOTAL_READS"

MIN_READS=$(echo "scale=2; $TOTAL_READS*$CUTOFF" | bc -l)
echo "$CUTOFF of total reads: $MIN_READS"

# Calculate read counts in negative control
THRESHOLD=$(grep ${NEG_CONTROL} read_counts.txt| awk '{print $2}')
echo "Number of reads for negative control ${NEG_CONTROL}: ${THRESHOLD}"

# filter read_counts.txt file for samples with more reads than negative control
# this code does not report samples which were filtered out
# awk -v limit="${THRESHOLD}" '$2 > limit' read_counts.txt 

# filter read_counts.txt file for samples with more reads than negative control
# AND reports samples which were filtered out
awk -v thresh="$THRESHOLD" -v col=2 '
{
    # Check if the value in the specified column is greater than the threshold
    # The +0 coerces the field value to a number to ensure a numerical comparison
    if ($col + 0 > thresh) {
        # Lines that meet the condition are redirected to a file
        print > "sample_read_count_above_nc.txt"
    } else {
        # Lines that do NOT meet the condition are printed to standard output
        print > "sample_read_count_below_nc.txt"
    }
}' "read_counts.txt"

NUM_EXCLUDED=$(cat sample_read_count_below_nc.txt | wc -l)
echo "Number of samples with read count lower than negative control: $NUM_EXCLUDED"
cat sample_read_count_below_nc.txt

# filter sample_read_count_above_nc.txt file for samples with more reads than specified cutoff
# AND reports samples which were filtered out
awk -v thresh="$MIN_READS" -v col=2 '
{
    # Check if the value in the specified column is greater than the threshold
    # The +0 coerces the field value to a number to ensure a numerical comparison
    if ($col + 0 > thresh) {
        # Lines that meet the condition are redirected to a file
        print > "sample_read_count_above_cutoff.txt"
    } else {
        # Lines that do NOT meet the condition are printed to standard output
        print > "sample_read_count_below_cutoff.txt"
    }
}' "sample_read_count_above_nc.txt"

NUM_EXCLUDED=$(cat sample_read_count_below_cutoff.txt | wc -l)
echo "Number of samples with read count lower than specified cutoff of $CUTOFF ($MIN_READS): $NUM_EXCLUDED"
cat sample_read_count_below_cutoff.txt

echo "Job finised with exit code $? at: `date`"
