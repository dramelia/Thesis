#!/bin/bash
#SBATCH --job-name=getGenomes
#SBATCH --time=04:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Run twice, once to get the fasta, and again to get the annotations.

# PATH = "my/file.path" - change this to the file path 

echo get EBV gtf

# link below  may be the wrong file, used old saved version from EBV_Matt folder in bioinformatics folders on my laptop and onedrive
# wget https://github.com/flemingtonlab/public/blob/master/annotation/chrEBV_Akata_inverted_refined_genes_plus_features_annotation_cleaned.gtf
# Use unzipped files for this


echo  cat_gtf start

cat PATH/gencode.v47.primary_assembly.annotation.gtf PATH/chrEBV_Akata_inverted_refined_genes_annotation_cleaned_GencodeL.gtf > PATH/gencode.v47.EBV1.annotation.gtf

echo cat_gtf finished
