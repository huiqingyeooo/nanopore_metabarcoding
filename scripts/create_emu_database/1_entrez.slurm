#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/entrez

TAX_ID=9895
TAXA_NAME=Bovidae
WORKING_FOLDER=/work/soghigian_lab/huiqing.yeo/metabarcoding/ref_database/coi

# Note: NCBI key is specific to user
# create your own using your NCBI account and paste it here
export NCBI_API_KEY=8d70027be2ed63ac8e322bdf53042857f608

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel -j 1 --link sbatch entrez.slurm ::: ${taxID} ::: ${taxaName} 

echo "starting run at: `date`"
####### script ##################################
cd ${WORKING_FOLDER}

# use just the esearch command to get an idea how many hits are found
# e.g., esearch -db nucleotide -query "txid7147[Organism:exp] AND (COI[Gene] OR COI[Gene]...)"

esearch -db nucleotide -query "txid${TAX_ID}[Organism:exp] \
AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene] \
OR cytochrome c oxidase subunit I[Gene] OR cytochrome c oxidase subunit 1[Gene] \
OR cytochrome oxidase subunit 1[Gene] OR cytochrome oxidase subunit I[Gene] \
OR COI[Title] OR CO1[Title] OR COX1[Title] OR COXI[Title] \
OR cytochrome c oxidase subunit I[Title] OR cytochrome c oxidase subunit 1[Title] \
OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title]) \
NOT complete genome[Title] \
AND 300:2000[SLEN]" \
| efetch -format uid  > uids_${TAX_ID}_${TAXA_NAME}.txt

#esearch -db nucleotide -query "txid${1}[Organism:exp] \
#AND (ITS2[Gene] OR its2 [Gene] \
#OR internal transcribed spacer 2[Gene] \
#OR internal transcribed spacers 1 and 2 [Gene] \
#OR 28S-18S intergenic spacer [Gene] \
#OR 28s-18s ribosomal RNA intergenic spacer [Gene] \
#OR 28s-18s rRNA intergenic spacer [Gene] \
#OR 28s-18s rRNA IGS[Gene] \
#OR ITS2[Title] OR its2 [Title] \
#OR internal transcribed spacer 2[Title] \
#OR internal transcribed spacers 1 and 2 [Title]
#OR 28S-18S intergenic spacer [Title] \
#OR 28s-18s ribosomal RNA intergenic spacer [Title] \
#OR 28s-18s rRNA intergenic spacer [Title] \
#OR 28s-18s rRNA IGS[Title]) \
#NOT complete genome[Title] \
#AND 300:2500[SLEN]" \
#| efetch -format uid  > uids_${1}_${2}.txt

echo "retrieval of uids for taxid ${TAX_ID} ${TAXA_NAME} done"

# note: there might be a curl command error
# check that number of sequences from esearch query and number of uids downloaded match
esearch -db nucleotide -query "txid${TAX_ID}[Organism:exp] \
AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene] \
OR cytochrome c oxidase subunit I[Gene] OR cytochrome c oxidase subunit 1[Gene] \
OR cytochrome oxidase subunit 1[Gene] OR cytochrome oxidase subunit I[Gene] \
OR COI[Title] OR CO1[Title] OR COX1[Title] OR COXI[Title] \
OR cytochrome c oxidase subunit I[Title] OR cytochrome c oxidase subunit 1[Title] \
OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title]) \
NOT complete genome[Title] \
AND 300:2000[SLEN]" | grep "Count" | sed "s#<Count>##g;s#</Count>##g;s#  ##g" > tmp_file

num_of_uids_detected=$(cat tmp_file)
num_of_uids_downloaded=$(cat ${WORKING_FOLDER}/uids_${TAX_ID}_${TAXA_NAME}.txt | wc -l)

# check that number matches

if [[ "$num_of_uids_detected" == "$num_of_uids_downloaded" ]]; then
    echo "Number of uids detected = Number of uids downloaded."
else
	echo "Number of uids detected not equals to Number of uids downloaded."
	exit 0
fi

# split uids text file to avoid 'too many requests' and 'connection reset' error
mkdir split_${TAXA_NAME}
cd split_${TAXA_NAME}

split -d -l 10000 ${WORKING_FOLDER}/uids_${TAX_ID}_${TAXA_NAME}.txt uids_${TAX_ID}_${TAXA_NAME}_chunk

echo "splitting of uids for taxid ${TAX_ID} ${TAXA_NAME} done"

for i in $(ls uids_${TAX_ID}_${TAXA_NAME}_chunk*);
do
chunkNo=`echo $i | awk -F"_" '{print $4}'`
cat ${i} | efetch -db nuccore -format fasta > sequences_${TAX_ID}_${TAXA_NAME}_${chunkNo}.fasta
echo "sequences_${TAX_ID}_${TAXA_NAME}_${chunkNo}.fasta download done"
done

echo "retrieval of sequences for taxid ${TAX_ID} ${TAXA_NAME} done"

cat sequences_${TAX_ID}_${TAXA_NAME}_*.fasta > sequences_${TAXA_NAME}_all.fasta

# note: there might be a curl command error
# check that number of sequences downloaded match counts

num_of_uids_downloaded=$(cat ${WORKING_FOLDER}/uids_${TAX_ID}_${TAXA_NAME}.txt | wc -l)
num_of_sequences_downloaded=$(grep ">" sequences_${TAXA_NAME}_all.fasta | wc -l)

# check that number matches
if [[ "$num_of_uids_downloaded" == "$num_of_sequences_downloaded" ]]; then
    echo "Number of uids downloaded = Number of sequences downloaded."
    cp sequences_${TAXA_NAME}_all.fasta ${WORKING_FOLDER}
    cd ${WORKING_FOLDER}
    rm -r split_${TAXA_NAME}
else
	echo "Number of uids downloaded not equals to Number of sequences downloaded."
	exit 0
fi

echo "Job finised with exit code $? at: `date`"
