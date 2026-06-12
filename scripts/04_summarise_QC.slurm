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

WORKING_DIR=${1}

# summarize nanoplot QC results
python run_QC_summary.py --input ${WORKING_DIR} --output ${WORKING_DIR}/summary

echo "Job finised with exit code $? at: `date`"
[huiqing.yeo@arc scripts]$ cat 04_summarise_QC.slurm 
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

WORKING_DIR=${1}

# summarize nanoplot QC results
python run_QC_summary.py --input ${WORKING_DIR} --output ${WORKING_DIR}/summary

echo "Job finised with exit code $? at: `date`"
