#!/bin/bash
#SBATCH --job-name=getGenomes
#SBATCH --time=04:00:00
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Run in folder where $1 (first argement on bottom line) is.

process_expression_table() {
    local input_file="$1"
    local output_filename="$2"
    local outdir="ForPGCNA"
    local output_file="${outdir}/${output_filename}"
    local fileinfo="${outdir}/#FileInfo.txt"

    echo $input_file

    # Make output directory if it doesn't exist
    mkdir -p "$outdir"

    # Process file
    awk 'BEGIN {OFS="\t"}
    NR==1 {
        printf "gene";
        for (i=4; i<=NF; i++) printf "%s%s", OFS, $i;
        printf "\n";
        next;
    }
    {
        gene = $1 "#" $3;
        printf "%s", gene;
        for (i=4; i<=NF; i++) printf "%s%s", OFS, $i;
        printf "\n";
    }' "$input_file" > "$output_file"

    # Write header to FileInfo.txt if it doesn't exist
    if [ ! -f "$fileinfo" ]; then
        echo -e "File\tHeaderLines" > "$fileinfo"
    fi

    # Add file info line
    echo -e "${output_filename}\t1" >> "$fileinfo"
}
# Process the file
process_expression_table "RNAseq_CLL_RSEM_Gene_count_VST_wGeneSymbols.tsv" "RNAseq_CLL_RSEM_Gene_count_VST_wGeneSymbols_PGCNA.txt"


