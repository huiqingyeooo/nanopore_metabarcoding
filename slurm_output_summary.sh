outputDir=/work/soghigian_lab/huiqing.yeo/metabarcoding/20260804_COI_MB005-MB016/0_fastq_files

# check slurm log files
for i in $(ls ${outputDir}/slurm-*.out);
do
col1=`echo ${i}`
#col2=`grep "sample_" ${i}`
col2=`grep -e "-pipiens" ${i}`
col3=`tail -n 1 ${i}`
echo $col1 $col2 $col3 >> ${outputDir}/slurm_summary.txt
done


