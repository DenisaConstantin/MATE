#!/bin/bash
# mate.sh - Main MATE Pipeline Script

GENOTYPE=""
REFERENCE=""
TE_SEQS=""
TE_COORD=""
OUT_DIR="mate_output"
#ALIAS_FILE=""


while [[ "$#" -gt 0 ]]; do
    case $1 in
        -genotype) GENOTYPE=$(realpath "$2"); shift ;;
        -reference) REFERENCE=$(realpath "$2"); shift ;;
        -te) TE_SEQS=$(realpath "$2"); shift ;;
        -te_coord) TE_COORD=$(realpath "$2"); shift ;;
        -out) OUT_DIR=$(realpath "$2"); shift ;;
#        -aliases) ALIAS_FILE=$(realpath "$2"); shift ;;
        *) echo "Eroare: Parametru necunoscut $1"; exit 1 ;;
    esac
    shift
done


if [ -z "$GENOTYPE" ] || [ -z "$REFERENCE" ] || [ -z "$TE_SEQS" ] || [ -z "$TE_COORD" ]; then
    echo "Usage: $0 -genotype <file.fasta> -reference <ref.fasta> -te <te.fasta> -te_coord <te_name_coord.txt> -out <output_folder>"
    exit 1
fi


export MATE_HOME=$(dirname $(realpath "$0"))


mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1


blastn -subject "$GENOTYPE" -query "$TE_SEQS" -outfmt 6 -out all_te_result.csv


awk -F'\t' '
{
    filename = $1
    gsub(/\//, "_", filename) 
    outfile = filename ".te.csv"
    print $0 >> outfile
}' all_te_result.csv


for i in *.te.csv; do 
    folder_name="${i%.te.csv}"
    mkdir -p "$folder_name"
    mv "$i" "$folder_name/"
done




for TE_DIR in */; do

    cd "$TE_DIR" || continue
    
  
    TE_NAME="${TE_DIR%/}"
 
    
   
    bash "$MATE_HOME/mate_core.sh" "$REFERENCE" "$GENOTYPE" "$TE_COORD"
    
    if [ -d "results" ]; then
        cd results || exit 1
        
        bash "$MATE_HOME/mate_unify.sh"
        bash "$MATE_HOME/mate_analyze.sh"
        
        
        cd ..
    else
        echo "Warning: The folder 'results' was not created in $TE_NAME!"
    fi
    
   
    cd ..
done

