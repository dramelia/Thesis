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


# THE DEFAULT FOR NF-RNASEQ IS SALMON!! USE _RSEM SCRIPT TO TEST RSEM 

# PATH = "my/file.path" - change this to a file path if required

nextflow run nf-core/rnaseq -r dev \
    --input ./SRA_files_Lopes_all/samplesheet/samplesheet_Lopes_edited.csv \
    --outdir ./RNAseq_files_Lopes_all \
    --gtf ./genomes/gencode.v48.EBV.primary_assembly.annotation.gtf \
    --fasta ./genomes/GRCh38.p14.EBV.genome.fa \
    --aligner star_rsem \
    --save_reference \
    -profile apptainer \
    -c nextflow.config_rnaseq \
    -resume