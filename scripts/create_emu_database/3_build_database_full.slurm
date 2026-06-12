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

database=/work/soghigian_lab/huiqing.yeo/filaroidDiv/emu/ref_database
refName=nem_cerv_hs_dip

echo "starting run at: `date`"
####### script ##################################

# concatenate coi sequences from desired taxa
cat sequences_Nematoda_all.fasta sequences_Cervidae_all.fasta \
sequences_Homo.sapiens_all.fasta sequences_Diptera_all.fasta \
> ${refName}.fasta

#sequences_Diptera_all.fasta 

# extract accession list
grep "^>" ${refName}.fasta \
| sed 's/^>//' \
| cut -d' ' -f1 \
> ${refName}.accessions

# concatenate taxids
cat \
sequences_Homo.sapiens.accessions.taxid \
sequences_Cervidae.accessions.taxid \
sequences_Nematoda.accessions.taxid \
sequences_Diptera.accessions.taxid \
> ${refName}.accessions.taxid

#sequences_Diptera.accessions.taxid

# sort taxid file to obtain the same order of accessions
TMP_FILE=$(mktemp)
while read LINE; do
    grep "$LINE" "${refName}.accessions.taxid" >>"$TMP_FILE"
done <"${refName}.accessions"

mv -f "$TMP_FILE" "${refName}_sorted.accessions"

### !!!!!
#this is much faster than the above code!
awk 'FNR == NR { order[$1] = NR; next } { print order[$1], $0 }' ${refName}.accessions ${refName}.accessions.taxid | sort -n | cut -d' ' -f2- > ${refName}_sorted.accessions
    
# download taxdump for nodes.dmp and names.dmp files
#conda deactivate
#conda activate /work/soghigian_lab/apps/conda/envs/taxonkit
#wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
#tar -zxvf taxdump.tar.gz -C taxdump
#conda deactivate
conda activate /work/soghigian_lab/apps/conda/envs/emu-3.6.2

# prepare seq2tax.map file
sed "s/ /_/g;s/\;//g;s/\,//g;s/\:/_/g;s#(##g;s#)##g;s/\./_/g;s/\-/_/g;s#[##g;s#]##g" ${refName}.fasta > ${refName}_format.fasta
grep ">" ${refName}_format.fasta > tmp1.txt
sed -i "s/>//g" tmp1.txt
cat ${refName}_sorted.accessions | awk -F"\t" '{print $2}' > tmp2.txt
paste -d'\t' tmp1.txt tmp2.txt > ${refName}_seq2tax.map
rm tmp*

# generate database
#emu build-database ref_nem_cerv --sequences ${refName}_format.fasta \
#--seq2tax ${refName}_seq2tax.map \
#--ncbi-taxonomy /work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database/taxdump

# remove lines from seq2tax and fasta files where taxid is not found in nodes.dmp
# temporary measure
# probably mismatch between ncbi taxid and taxonomy databases for scutellonema afribrachyurus
grep -v "3686840" ${refName}_seq2tax.map > ${refName}_seq2tax_2.map
grep "3686840" ${refName}_seq2tax.map | awk -F"\t" '{print $1}' > tmp.txt
grep ">" ${refName}_format.fasta | sed "s/>//g" | grep -vxf tmp.txt > tmp_subset.txt
seqtk subseq ${refName}_format.fasta tmp_subset.txt > ${refName}_format_2.fasta

# generate database
/work/soghigian_lab/apps/conda/envs/emu-3.6.2/emu build-database ${refName} --sequences ${refName}_format_2.fasta \
--seq2tax ${refName}_seq2tax_2.map \
--ncbi-taxonomy ${database}/taxdump

echo "Job finised with exit code $? at: `date`"
