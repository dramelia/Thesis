#!/bin/bash
#SBATCH --job-name=testscript1	# Job name
#SBATCH --time=01:00:00		# Request runtime (hh:mm:ss)
#SBATCH --mem=25G		# Request memory
#SBATCH --ntasks=1		# Number of tasks
#SBATCH --cpus-per-task=8	# Number of CPUs. If I need GPUs see Aire docs for code.

# Load ny necessary modules
# module load <module name>

# Execute your application
# ./example.bin		this runs an example binary file in the current directory.
# Put my script here and it will run in the current directory.

# All of abiove will be ignored if running via nextflow

# This script uses the default nextflow.config file

# PATH = "my/file.path" - change this to the file path 

nextflow run nf-core/fetchngs -r dev \
--input PATH/GSE199868_Lopes_SRR_Acc_List.csv \
--outdir PATH/SRA_files_Lopes_all \
--nf_core_pipeline rnaseq \
-c nextflow.config_sraNoSlurm