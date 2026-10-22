#!/bin/bash
# mate_associate.sh

INPUT_INS=$1
COORD_TXT=$2

awk '

NR==FNR {
    chrom[NR] = $1
    start[NR] = $2
    end[NR] = $3
    name[NR] = $4
    max_records = NR
    next
}

{
    id = $1
    chrom1 = $2
    coord1 = $16
    coord2 = $17

    if (coord1 !~ /^[0-9]+$/) next

  
    for (i = 1; i <= max_records; i++) {
        if (chrom1 == chrom[i]) {
            if ((coord1 >= start[i] && coord1 <= end[i]) || (coord2 >= start[i] && coord2 <= end[i])) {
                print id "\t" name[i]
            }
        }
    }
}
' "$COORD_TXT" "$INPUT_INS"
