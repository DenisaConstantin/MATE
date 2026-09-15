#!/bin/bash
# mate.sh - Main MATE Pipeline Script

# Initializare variabile
GENOTYPE=""
REFERENCE=""
TE_SEQS=""
TE_COORD=""
OUT_DIR="mate_output"
#ALIAS_FILE=""

# Parsare argumente
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

# Verificare argumente obligatorii
if [ -z "$GENOTYPE" ] || [ -z "$REFERENCE" ] || [ -z "$TE_SEQS" ] || [ -z "$TE_COORD" ]; then
    echo "Usage: $0 -genotype <file.fasta> -reference <ref.fasta> -te <te.fasta> -te_coord <te_name_coord.txt> -out <output_folder>"
    exit 1
fi

# Salvarea caii catre folderul unde se afla toate scripturile MATE
#MATE_HOME=$(dirname $(realpath "$0"))
export MATE_HOME=$(dirname $(realpath "$0"))

# Creare si mutare in folderul de output
mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

#echo "[MATE] 1/5: Rulare blastn genotip vs transpozoni..."
blastn -subject "$GENOTYPE" -query "$TE_SEQS" -outfmt 6 -out all_te_result.csv

#echo "[MATE] 2/5: Impartire rezultate per transpozon..."
awk -F'\t' '
{
    filename = $1
    gsub(/\//, "_", filename) 
    outfile = filename ".te.csv"
    print $0 >> outfile
}' all_te_result.csv

#echo "[MATE] Creare foldere specifice pentru transpozoni..."
for i in *.te.csv; do 
    folder_name="${i%.te.csv}"
    mkdir -p "$folder_name"
    mv "$i" "$folder_name/"
done

#echo "[MATE] Incepem procesarea individuala pentru fiecare transpozon..."

# 3, 4 & 5: Iteram prin fiecare folder creat (fiecare transpozon)
for TE_DIR in */; do
    # Intram in folderul transpozonului (ex: hobo/)
    cd "$TE_DIR" || continue
    
    # Inlaturam slash-ul de la final pentru afisare curata
    TE_NAME="${TE_DIR%/}"
    #echo "=================================================="
    #echo "[MATE] Procesam transpozonul: $TE_NAME"
    #echo "=================================================="
    
    #echo "[MATE] 3/5: Rulare mapare principala (mate_core.sh)..."
    bash "$MATE_HOME/mate_core.sh" "$REFERENCE" "$GENOTYPE" "$TE_COORD"
    
    #echo "[MATE] 4/5 si 5/5: Rulare unificare si analiza..."
    if [ -d "results" ]; then
        # Intram in folderul results specific acestui transpozon
        cd results || exit 1
        
        bash "$MATE_HOME/mate_unify.sh"
        bash "$MATE_HOME/mate_analyze.sh"
        
        # Iesim din results inapoi in folderul transpozonului
        cd ..
    else
        echo "Avertisment: Folderul 'results' nu a fost creat in $TE_NAME!"
    fi
    
    # Iesim din folderul transpozonului inapoi in OUT_DIR pentru a continua bucla
    cd ..
done

#echo "=================================================="
#echo "[MATE] Analiza s-a incheiat cu succes pentru toti transpozonii!"
#echo "Rezultatele finale se afla in: $OUT_DIR/<nume_transpozon>/results"
