#!/bin/bash
#SBATCH --job-name=getGenomes
#SBATCH --time=06:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# EBV version
# This contains my edits of the FlemingtonLabGIT files (from https://github.com/flemingtonlab/public). 

## Background
# Want to combine the Human Gencode files with the EBV files from the FlemingtonLabGIT repository.
    
## Files:

# chrEBV_Akata_inverted_2.fa : Original file from FlemingtonLabGIT/annotation
# chrEBV_Akata_inverted_refined_genes_annotation_cleaned_GencodeL.gtf : converted to gencode-like format
# GRCh38.p14.genome.fa.gz
# gencode.v48.annotation.gtf.gz
# gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz
# gencode.v48.primary_assembly.annotation.gtf.gz

## Commands run

# GTF file 1*

#cat gencode.v48.annotation.gtf chrEBV_Akata_inverted_refined_genes_annotation_cleaned_GencodeL.gtf > gencode.v48.EBV.annotation.gtf

# GTF file 2*

#cat gencode.v48.primary_assembly.annotation.gtf chrEBV_Akata_inverted_refined_genes_annotation_cleaned_GencodeL.gtf > gencode.v48.EBV.primary_assembly.annotation.gtf

# GTF file 3*

#cat gencode.v48.chr_patch_hapl_scaff.annotation.gtf chrEBV_Akata_inverted_refined_genes_annotation_cleaned_GencodeL.gtf > gencode.v48.EBV.chr_patch_hapl_scaff.annotation.gtf

# FA file*

cat GRCh38.p14.genome.fa chrEBV_Akata_inverted_2.fa > GRCh38.p14.EBV.genome.fa


