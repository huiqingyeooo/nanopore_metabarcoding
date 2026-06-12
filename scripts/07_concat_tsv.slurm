#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=50G

####### Set environment variables ###############
PID=${1}
LEN=${2}
PARENT_DIR=${3}/pid${PID}_len${LEN}
DEMULTIPLEX_DIR=${4}
SAMPLE_INFO=${5}

echo "starting run at: `date`"
####### script ##################################
cd ${PARENT_DIR}

# copies out rel-abundance.tsv from each sample folder

for dir in "$PARENT_DIR"/*/; do
    # Remove trailing slash and extract subfolder name
    subfolder_name=$(basename "$dir")
    source_file="$dir/${subfolder_name}_trimmed.fastq_rel-abundance.tsv"
    dest_file="$PARENT_DIR/${subfolder_name}_trimmed.fastq_rel-abundance.tsv"

    # Check if results.csv exists in the subdirectory
    if [[ -f "$source_file" ]]; then
        cp "$source_file" "$dest_file"
        echo "Copied $source_file to $dest_file"
    else
        echo "No results.csv in $dir"
    fi
done

########################################
# generate species-abundance matrix
########################################
# wrangle dataframe for generating species-abundance matrix
for file in *.tsv; do
  sed -i "s/ /_/g" $file
  sample=$(basename $file _trimmed.fastq_rel-abundance.tsv)
  echo "processing ${sample}"
  
  # print only relevant columns; taxid, species names, and abundance
  awk -F"\t" '
  NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
        }
    }
    { print $(f["tax_id"])"_"$(f["species"])"\t"$(f["estimated_counts"]) }' $file > $file.tmp
    
    sed -i '1d' $file.tmp

  # transpore dataframe
  awk '{
  for (i=1; i<=NF; i++) {
    a[i] = (a[i] ? a[i] FS $i : $i)
    }}
    END {
    for (i in a) {
    print a[i]
    }}' $file.tmp > $file.tmp2
    
  # Add header 'sample' into first columm, first row
  # Add sample name ($sample) into first column, second row
  # Search and replace ' ' with '\t'
    sed '1s#^#sample\t#' $file.tmp2 | sed "2s#^#$sample\t#" | sed "s# #\t#g"> $file.tmp3
done

# generate species-abundance matrix
source ~/software/init-conda-sl #activate python 3

LOC_PYTHON_SCRIPT=/work/soghigian_lab/huiqing.yeo/metabarcoding/scripts

python ${LOC_PYTHON_SCRIPT}/generate_species_abundance_matrix.py \
"*_trimmed.fastq_rel-abundance.tsv.tmp3" \
species_abundance_matrix.csv


########################################
# summarise top two hits per sample
########################################

for file in *.tsv; do
sample=$(basename $file _trimmed.fastq_rel-abundance.tsv)
TOTAL_READS=`awk -F"\t" '{sum += $2} END {print sum}' ${file}.tmp`
sort -nk 2,2 -nr ${file}.tmp > ${file}.tmp4
TOP_1_HIT=`sed -n '1p' ${file}.tmp4`
TOP_2_HIT=`sed -n '2p' ${file}.tmp4`
TOP_3_HIT=`sed -n '3p' ${file}.tmp4`
echo -e "$sample\t$TOTAL_READS\t$TOP_1_HIT\t$TOP_2_HIT\t$TOP_3_HIT"  >> top_hits.tmp
done

# calculate ratio of top1/top2 hit
awk -F"\t" '{print $0, $4/$6}' OFS="\t" top_hits.tmp | \
awk -F"\t" '{if ($9 >= 5) $10 = "successful"; else $10 = "unsuccessful"; print}' OFS="\t" > top_hits.csv

# add header into long form csvcat
header="sample\ttotal_reads\ttop_1_ID\ttop_1_count\ttop_2_ID\ttop_2_count\ttop_3_ID\ttop_3_count\ttop1/top2_ratio\tstatus"
sed -i "1i\\$header" top_hits.csv

####################################################
# concatenate rel-abundance.tsvs into long form csv
####################################################
for file in *.tsv; do
sample=$(basename $file _trimmed.fastq_rel-abundance.tsv)
sed '1d' $file | \
awk -v var="$sample" '{print var, $0}' OFS="\t" > $file.tmp5
done

cat *.tmp5 > rel_abundance.csv

# add header into long form csvcat
header="sample\ttax_id\tabundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\tsubspecies\tspecies_subgroup\tspecies_group\testimated_counts"
sed -i "1i\\$header" rel_abundance.csv

rm *tmp*
rm *tsv*

######################################
# match results with sample_info.csv
######################################
TOP_HITS_HEADER="total_reads\ttop_1_ID\ttop_1_count\ttop_2_ID\ttop_2_count\ttop_3_ID\ttop_3_count\ttop1/top2_ratio\tstatus"
SAMPLE_INFO_HEADER=`cat ${DEMULTIPLEX_DIR}/${SAMPLE_INFO} | head -n 1`
HEADER="$SAMPLE_INFO_HEADER\t$TOP_HITS_HEADER"

sed '1d' ${DEMULTIPLEX_DIR}/${SAMPLE_INFO} > ${SAMPLE_INFO}.tmp

awk '{ split($1, a, "_"); $1=a[1]; print }' OFS="\t" top_hits.csv | \
sed '1d' > top_hits.csv.tmp

# left join, first file is the full list, second file is results
join -a 1 -t $'\t' -e 'NULL' -o auto <(sort ${SAMPLE_INFO}.tmp) <(sort top_hits.csv.tmp) > results.csv
sed -i "1i\\$HEADER" results.csv
rm *tmp*

echo "Job finised with exit code $? at: `date`"
