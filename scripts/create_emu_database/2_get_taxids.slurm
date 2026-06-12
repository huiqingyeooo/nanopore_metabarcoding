#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/entrez

PATH=/work/soghigian_lab/apps/seqtk/:$PATH
database=/work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 sbatch 2_get_taxids.slurm ::: Cervidae

echo "starting run at: `date`"
####### script ##################################
cd ${database}

echo $SLURM_JOBID

# get taxids for each taxa ref database
# extract accession numbers
#grep "^>" sequences_${1}_all.fasta \
#| sed 's/^>//' \
#| cut -d' ' -f1 \
#> sequences_${1}.accessions

# get taxids from ncbi using list of accession numbers
#split -d -l 1000 sequences_${1}.accessions chunk_${1}
#chunk_files=$(ls chunk_${1}*)

export NCBI_API_KEY=8d70027be2ed63ac8e322bdf53042857f608

#for i in ${chunk_files};
#do
#echo getting taxids from ${i}
#epost -db nucleotide -input ${i} \
#| esummary \
#| xtract -pattern DocumentSummary -element AccessionVersion TaxId Organism > sequences_${1}.accessions.${i}.taxid
#done

#cat sequences_${1}.accessions.chunk_${1}*.taxid > sequences_${1}.accessions.taxid

# check that number of taxids and sequences matches
#num_of_taxids=$(cat sequences_${1}.accessions.taxid | wc -l)
#num_of_sequences=$(grep ">" sequences_${1}_all.fasta | wc -l)

#if [[ "$num_of_taxids" == "$num_of_sequences" ]]; then
#    echo "Number of taxids = Number of sequences."
#    rm *chunk*
#else
#	echo "Number of taxids not equals to Number of sequences."
#	rm _chunk_summary.txt
#	for i in chunk_*;
#	do
#	col1=`cat ${i} | wc -l`
#	col2=`cat sequences_${1}.accessions.${i}.taxid | wc -l`
#	echo ${i} $col1 $col2 >> _chunk_summary.txt
#	done
#	awk '$2 != $3 {print $1}' _chunk_summary.txt > _chunk_to_rerun.txt
#fi	


# Troubleshoot entrez errors when retrieving taxid

# combine taxids
for i in $(cat _chunk_to_rerun.txt)
do
echo getting taxids from ${i}
epost -db nucleotide -input ${i} \
| esummary \
| xtract -pattern DocumentSummary -element AccessionVersion TaxId Organism > sequences_${1}.accessions.${i}.taxid
done

rm sequences_${1}.accessions.taxid
cat sequences_${1}.accessions.chunk_${1}*.taxid > sequences_${1}.accessions.taxid

# check that number of taxids and sequences matches
num_of_taxids=$(cat sequences_${1}.accessions.taxid | wc -l)
num_of_sequences=$(grep ">" sequences_${1}_all.fasta | wc -l)

if [[ "$num_of_taxids" == "$num_of_sequences" ]]; then
    echo "Number of taxids = Number of sequences."
    #rm *chunk*
else
	echo "Number of taxids not equals to Number of sequences."
	rm _chunk_summary.txt
	for i in chunk_*;
	do
	col1=`cat ${i} | wc -l`
	col2=`cat sequences_Diptera.accessions.${i}.taxid | wc -l`
	echo ${i} $col1 $col2 >> _chunk_summary.txt
	done
	awk '$2 != $3 {print $1}' _chunk_summary.txt > _chunk_to_rerun.txt
fi	

# combine taxids
#for i in $(cat _chunk_to_rerun.txt)
#do
#echo getting taxids from ${i}
#epost -db nucleotide -input ${i} \
#| esummary \
#| xtract -pattern DocumentSummary -element AccessionVersion TaxId Organism > sequences_${1}.accessions.${i}.taxid
#done

echo "Job finised with exit code $? at: `date`"
