#!/bin/bash
# mate_associate.sh - Asociere directa pe baza de coordonate (optimizat in AWK)

INPUT_INS=$1
COORD_TXT=$2

awk '
# 1. Citim mai întâi fișierul de referință (COORD_TXT) și îi reținem datele
NR==FNR {
    chrom[NR] = $1
    start[NR] = $2
    end[NR] = $3
    name[NR] = $4
    max_records = NR
    next
}
# 2. Citim fișierul de inserții și aplicăm exact aceleași verificări
{
    id = $1
    chrom1 = $2
    coord1 = $16
    coord2 = $17

    # Verificare numerică (echivalentul lui: [[ "$coord1" =~ ^[0-9]+$ ]] || continue)
    if (coord1 !~ /^[0-9]+$/) next

    # Căutare în setul de referință stocat
    for (i = 1; i <= max_records; i++) {
        if (chrom1 == chrom[i]) {
            # Verificarea intervalelor pentru coord1 sau coord2
            if ((coord1 >= start[i] && coord1 <= end[i]) || (coord2 >= start[i] && coord2 <= end[i])) {
                print id "\t" name[i]
            }
        }
    }
}
' "$COORD_TXT" "$INPUT_INS"
