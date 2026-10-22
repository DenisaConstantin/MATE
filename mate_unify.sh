#!/bin/bash

INPUT="unannotated_conserved_ins.csv"
OUTPUT="unified_unannotated_conserved_insertions.tsv"
THRESHOLD=50  


sort -k1,1 -k2,2 -k14,14n "$INPUT" | awk -v thresh="$THRESHOLD" '
function abs(x) { return (x < 0) ? -x : x }


function print_record() {
    if (saved_line == "") return;
    #Splituim linia originala pentru a-i actualiza doar coordonatele
    split(saved_line, fields, " ");
    fields[14]=c_qs; fields[15]=c_qe; fields[16]=c_rs; fields[17]=c_re;
    
    for(i=1; i<=17; i++) printf "%s%s", fields[i], (i==17 ? "" : "\t");
    print "";
}

BEGIN { FS=" "; OFS="\t"; saved_line=""; }

{
    
    split($1, name_parts, "_");
    scaffold_base = name_parts[1];

    q_s=$14; q_e=$15; r_s=$16; r_e=$17; r_chr=$2;
    strand = (r_s < r_e) ? 1 : -1;

    if (saved_line == "") {
        
        c_scaff=scaffold_base; c_rchr=r_chr; c_qs=q_s; c_qe=q_e; c_rs=r_s; c_re=r_e; c_st=strand;
        saved_line=$0;
        next;
    }

  
    gap_ass = abs(q_s - c_qe);

  
    if (c_st == 1) gap_ref = abs(r_s - c_re);
    else gap_ref = abs(r_e - c_rs);

 
    if (scaffold_base == c_scaff && r_chr == c_rchr && strand == c_st && gap_ass <= thresh && gap_ref <= thresh) {
        # Actualizam coordonatele de final (unim fragmentul)
        if (q_e > c_qe) c_qe = q_e;
        if (c_st == 1) {
            if (r_e > c_re) c_re = r_e;
        } else {
            # Pe catena minus, r_s creste in sens invers
            if (r_s > c_rs) c_rs = r_s;
            if (r_e < c_re) c_re = r_e;
        }
    } else {
      
        print_record();
        c_scaff=scaffold_base; c_rchr=r_chr; c_qs=q_s; c_qe=q_e; c_rs=r_s; c_re=r_e; c_st=strand;
        saved_line=$0;
    }
}
END { print_record(); }
' > "$OUTPUT"

