import pandas as pd
import glob
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Concatenate files using glob pattern")
parser.add_argument("pattern", help="Glob pattern for input files")
parser.add_argument("output", help="Output file name")

# Parse arguments from command line
args=parser.parse_args()

# Get all matching files
#files = glob.glob("*_trimmed.fastq_rel-abundance.tsv.tmp3")
files = sorted(glob.glob(args.pattern))

if not files:
	raise ValueError("No files matched the given pattern")

# Read and store DataFrames in a list
dfs = [pd.read_csv(file, sep="\t") for file in files]
# Concatenate all DataFrames
abundance_matrix = pd.concat(dfs, axis=0, join='outer').fillna(0)

# Save
#abundance_matrix.to_csv("species_abundance_matrix.tsv", sep="\t", index=False)
abundance_matrix.to_csv(args.output, sep="\t", index=False)

print(f"Combined {len(files)} files into {args.output}")

# code for concatenating individual files (for troubleshooting)
#import pandas as pd
#data1=pd.read_csv('sample01_q18_trimmed.fastq_rel-abundance.tsv.tmp3',sep="\t")
#data2=pd.read_csv('sample02_q18_trimmed.fastq_rel-abundance.tsv.tmp3',sep="\t")
#data3=pd.read_csv('sample03_q18_trimmed.fastq_rel-abundance.tsv.tmp3',sep="\t")
#data4=pd.read_csv('sample06_q18_trimmed.fastq_rel-abundance.tsv.tmp3',sep="\t")
#abundance_matrix = pd.concat([data1,data2,data3,data4], axis=0, join='outer').fillna(0)
#print(abundance_matrix)
#abundance_matrix.to_csv("species_abundance_matrix.tsv", sep="\t", index=False)
