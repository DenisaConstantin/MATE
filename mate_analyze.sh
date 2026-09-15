#!/bin/bash

INPUT="unified_unannotated_conserved_insertions.tsv"
OUTPUT="unannotated_conserved_insertions.csv"
MARGIN=20
MIN_LEN=300  # Pragul minim de 300 nt pentru lungimea insertiei

#Header pentru fisierul final
echo -e "Coloane_Originale...\tDistanta_Scaffold\tTip_Insertie" > "$OUTPUT"

#Sortare: 1.Cromozom, 2.Start_Ref, 3.Scaffold, 4.Start_Query
sort -k2,2 -k16,16n -k1,1 -k14,14n "$INPUT" | awk -v m="$MARGIN" -v min_l="$MIN_LEN" '
function abs(x) { return (x < 0) ? -x : x }

#Extrage numele scaffold-ului (ex: din tig00001889:+)
function get_scaffold(id_str) {
    sub(/:[+-]$/, "", id_str)
    sub(/_[0-9]+-[0-9]+$/, "", id_str)
    return id_str
}

#Extrage catena (+ sau -)
function get_strand(id_str) {
    n = split(id_str, a, ":")
    return a[n]
}

BEGIN { 
    FS="\t"; OFS="\t"
    count = 0
}

{
    #Logica de grupare (daca e aceeasi regiune in referinta)
    if (NR == 1 || $2 != last_chrom || abs($16 - last_ref_start) > m || abs($17 - last_ref_end) > m) {
        if (NR > 1) process_group()
        
        count = 0
        last_chrom = $2
        last_ref_start = $16
        last_ref_end = $17
    }

    count++
    full_lines[count] = $0
    ids[count] = $1
    q_starts[count] = $14
    q_ends[count] = $15
}

END {
    if (count > 0) process_group()
}

function process_group() {
    #Filtrare: calculam lungimea pe referinta
    #len = abs(last_ref_end - last_ref_start)
     len = abs(q_ends[1] - q_starts[1])
    if (len < min_l) return;

    if (count == 1) {
        print full_lines[1], "N/A", "Unique_Insertion_Resolved"
        return
    }

    #1. Rezolvarea redundantei: alegem scaffold-ul cu aliniamentul cel mai lung
    best_idx = 1
    max_len = abs(q_ends[1] - q_starts[1])
    
    for (i = 2; i <= count; i++) {
        curr_len = abs(q_ends[i] - q_starts[i])
        if (curr_len > max_len) {
            max_len = curr_len
            best_idx = i
        }
    }
    best_scaff = get_scaffold(ids[best_idx])

    #2. Filtram doar insertiile care apartin scaffold-ului "castigator"
    new_count = 0
    for (i = 1; i <= count; i++) {
        if (get_scaffold(ids[i]) == best_scaff) {
            new_count++
            f_lines[new_count] = full_lines[i]
            f_ids[new_count] = ids[i]
            f_qs[new_count] = q_starts[i]
            f_qe[new_count] = q_ends[i]
        }
    }

    #3. Clasificarea finala a grupului curatat
    if (new_count == 1) {
        print f_lines[1], "N/A", "Unique_Insertion_Resolved"
        return
    }

    for (i = 1; i <= new_count; i++) {
        curr_strand = get_strand(f_ids[i])
        tag = ""
        dist = "N/A"

        if (i == 1) {
            tag = "First_in_Cluster"
            dist = 0
        } else {
            prev_strand = get_strand(f_ids[i-1])
            if (curr_strand == prev_strand) {
                dist = abs(f_qs[i] - f_qe[i-1])
                if (dist < 2000) tag = "Tandem_Cluster"
                else tag = "Dispersed_Copies"
            } else {
                tag = "Inverted_On_Same_Scaff"
            }
        }
        print f_lines[i], dist, tag
    }
}
' >> "$OUTPUT"
rm unannotated_conserved_ins.csv unified_unannotated_conserved_insertions.tsv
sed -i 's/\t/,/g' unannotated_conserved_insertions.csv
sed -i '1d' unannotated_conserved_insertions.csv
