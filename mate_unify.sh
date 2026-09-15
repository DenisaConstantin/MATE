#!/bin/bash

INPUT="unannotated_conserved_ins.csv"
OUTPUT="unified_unannotated_conserved_insertions.tsv"
THRESHOLD=50  #Distanta maxima (nucleotide) intre fragmente

#Procesare:
#Sortam prima data datele pentru a aduce fragmentele vecine unul langa altul
#Sorteaza dupa: Tig/Scaffold (Col 1), Ref Chr (Col 2), Start Asamblare (Col 14)
sort -k1,1 -k2,2 -k14,14n "$INPUT" | awk -v thresh="$THRESHOLD" '
function abs(x) { return (x < 0) ? -x : x }

#Functie care printeaza înregistrarea (fie ea unificata sau singura)
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
    #Extragem numele de baza al scaffold-ului (ex: din tig00001889_... luăm tig00001889)
    #Cautam prefixul pana la primul underscore sau folosim tot ID-ul daca nu are
    split($1, name_parts, "_");
    scaffold_base = name_parts[1];

    q_s=$14; q_e=$15; r_s=$16; r_e=$17; r_chr=$2;
    strand = (r_s < r_e) ? 1 : -1;

    if (saved_line == "") {
        # Prima linie din fișier
        c_scaff=scaffold_base; c_rchr=r_chr; c_qs=q_s; c_qe=q_e; c_rs=r_s; c_re=r_e; c_st=strand;
        saved_line=$0;
        next;
    }

    #Calculam distanta in asamblare (Assembly Gap)
    gap_ass = abs(q_s - c_qe);

    #Calculam distanta in referinta (Reference Gap)
    if (c_st == 1) gap_ref = abs(r_s - c_re);
    else gap_ref = abs(r_e - c_rs);

    #CONDITIA DE UNIRE (foarte specifica):
    #1. Acelasi scaffold de baza (tigXXXX)
    #2. Acelasi cromozom de referinta (2R, 3L etc)
    #3. Aceeasi orientare (plus/minus)
    #4. Gap-ul in ambele este sub pragul ales (50 nt)
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
        # Nu se pot unifica -> printam ce am salvat anterior si incepem un grup nou
        print_record();
        c_scaff=scaffold_base; c_rchr=r_chr; c_qs=q_s; c_qe=q_e; c_rs=r_s; c_re=r_e; c_st=strand;
        saved_line=$0;
    }
}
END { print_record(); }
' > "$OUTPUT"

