#!/bin/bash
#SBATCH --job-name=getGenomes
#SBATCH --time=04:00:00
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Run twice, once to get the fasta, and again to get the annotations.

echo getting EBV genome

# wget https://github.com/flemingtonlab/public/blob/master/annotation/chrEBV_Akata_inverted_2.fa

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Now swap over to getting the annotation

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz

# PATH = "my/file.path" - change this to the file path 

echo start cat_fa

cat PATH/GRCh38.primary_assembly.genome.fa /users/medafisa/scratch/PGCNA_work/genomes/chrEBV_Akata_inverted_2.fa > PATH/GRCh38.primary_assembly_EBV1.genome.fa

echo finished cat_fa

