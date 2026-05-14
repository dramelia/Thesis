#!/bin/bash
#SBATCH --job-name=getGenomes
#SBATCH --time=04:00:00
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Run twice, once to get the fasta, and again to get the annotations.

echo getting genome

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.p14.genome.fa.gz

# Now swap over to getting the annotation

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz

# Get the patch

`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz`

# Get the primary assembly annotation

`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz`


