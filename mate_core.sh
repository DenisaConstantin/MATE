
awk -F "\t" '{print > $2".ctg.csv"}' *.te.csv

DIST=50

for file in *.ctg.csv; do
    awk -v DIST="$DIST" -F'\t' '
    {
        line[NR] = $0
        sstart[NR] = $9 + 0
        send[NR]   = $10 + 0
        score[NR]  = $12 + 0
        keep[NR] = 1
    }

    function abs(x) { return x < 0 ? -x : x }
    function min(x, y) { return x < y ? x : y }
    function max(x, y) { return x > y ? x : y }

    END {
        for (i = 1; i <= NR; i++) {
            for (j = 1; j <= NR; j++) {
                if (i == j) continue

                smin_i = min(sstart[i], send[i])
                smax_i = max(sstart[i], send[i])
                smin_j = min(sstart[j], send[j])
                smax_j = max(sstart[j], send[j])

                if (smax_i == smin_j && abs(smin_i - smax_j) <= DIST) {
                    if (score[i] < score[j]) keep[i] = 0
                }

                if (smin_i == smax_j && abs(smax_i - smin_j) <= DIST) {
                    if (score[i] < score[j]) keep[i] = 0
                }

                if (smin_i >= smin_j && smax_i <= smax_j && score[i] < score[j]) {
                    keep[i] = 0
                }
            }
        }

        for (i = 1; i <= NR; i++) if (keep[i]) print line[i]
    }
    ' "$file" > "${file%.csv}_filtered.csv"
done

for f in *_filtered.csv; do
  awk '
  BEGIN { OFS="\t" }
  {
    line[NR] = $0
    query[NR] = $1
    contig[NR] = $2
    qstart[NR] = $7
    qend[NR] = $8
    sstart[NR] = $9
    send[NR] = $10
    bitscore[NR] = $12
  }
  END {
    for (i = 1; i <= NR; i++) {
      if (i in processed) continue

      matched = 0
      for (j = 1; j <= NR; j++) {
        if (i == j || (j in processed)) continue

        if (query[i] == query[j] && contig[i] == contig[j] &&
            send[i] == sstart[j]) {

          if ((sstart[i] < send[i]) && (sstart[j] < send[j])) {
            # sens crescător
            min_q = (qstart[i] < qstart[j]) ? qstart[i] : qstart[j]
            max_q = (qend[i] > qend[j]) ? qend[i] : qend[j]
            min_s = (sstart[i] < sstart[j]) ? sstart[i] : sstart[j]
            max_s = (send[i] > send[j]) ? send[i] : send[j]
          } else if ((sstart[i] > send[i]) && (sstart[j] > send[j])) {
          
            min_q = (qstart[i] < qstart[j]) ? qstart[i] : qstart[j]
            max_q = (qend[i] > qend[j]) ? qend[i] : qend[j]
            max_s = (sstart[i] > sstart[j]) ? sstart[i] : sstart[j]
            min_s = (send[i] < send[j]) ? send[i] : send[j]
          } else {
            continue
          }

          split(line[i], fields, "\t")
          fields[7] = min_q
          fields[8] = max_q
          fields[9] = min_s
          fields[10] = max_s
          print fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9], fields[10], fields[11], fields[12]

          processed[i] = 1
          processed[j] = 1
          matched = 1
          break
        }
      }
      if (!matched) print line[i]
    }
  }
  ' "$f" > "${f%.ctg.csv}_merged.csv"
done
rm *.ctg.csv *_filtered.csv

for i in *_merged.csv; do
  > "$i.ins.csv"
  while IFS=$'\t' read -r te tig _ _ _ _ qstart qend sstart send _ _; do
        ins=$((qend - qstart))
    if [ "$sstart" -lt "$send" ]; then
      strand="+"
      echo -e "$tig\t$sstart\t$send\t$strand\t$ins" >> "$i.ins.csv"
    elif [ "$sstart" -gt "$send" ]; then
      strand="-"
      echo -e "$tig\t$send\t$sstart\t$strand\t$ins" >> "$i.ins.csv"
    else
      strand="?"
    fi          
                
    
  done < "$i"                                                    
done

rm *_merged.csv
cat *.ins.csv > output_ins.txt
rm *.ins.csv

#Extrag subsecvente
while read tig start end strand ins
do 
#Catena plus
if [ $strand = + ] 
then 
startfl=1
if ! [ $start -le 3000 ]
then
startfl=$( expr $start - 3001 ) 
fi
endfl=$( expr $end + 3000 )
        echo -e "$tig\t$startfl\t$endfl\t$ins\t$start\t"+"\t$end" >> plus.bed

elif [ $strand = - ]
then 
startfl=1
if ! [ $start -le 3000 ]
then
startfl=$( expr $start - 3001 )
fi
endfl=$( expr $end + 3000 )
        echo -e "$tig\t$endfl\t$startfl\t$ins\t$end\t"-"\t$start" >> minus.bed
fi
done < output_ins.txt


seqkit subseq --bed plus.bed $2 > secv_plus.fasta
seqkit subseq --bed minus.bed $2 > secv_minus.fasta

blastn -subject $1 -query secv_plus.fasta -outfmt 7 -out plus_blast.csv
blastn -subject $1 -query secv_minus.fasta -outfmt 7 -out minus_blast.csv

for i in *_blast.csv
do
    output_dir="output_${i%_blast.csv}"
    mkdir -p "$output_dir"
    
  
    while IFS= read -r line; do
        if [[ $line == \#\ Query:* ]]; then
            #query_name=$(echo "$line" | sed 's/# Query: //; s/[^a-zA-Z0-9._-]/_/g')
            query_name=$(echo "$line" | sed 's/# Query: //; s/:.*//; s/[^a-zA-Z0-9._-]/_/g')
            current_file="${output_dir}/${query_name}.txt"
            echo "$line" > "$current_file"
        elif [[ -n $current_file ]]; then
            echo "$line" >> "$current_file"
        fi
    done < "$i"
done

cat plus.bed minus.bed > all.csv
sed -i 's/\t/ /g' all.csv
awk '{
    if ($6 == "-") {
        tmp = $5; $5 = $7; $7 = tmp;
        tmp = $2; $2 = $3; $3 = tmp;
    }
    print
}' all.csv > all_ins.csv
rm all.csv

current_dir=$(pwd)

for dirname in ./*/
do
  cd "$dirname" || continue


for f in *.txt; do echo ${f%.txt} >> insertii; sed -i 's/\t/ /g' $f; awk 'NR > 4 {print $0}' $f > ${f%.txt}".int.txt"; done


for i in *.int.txt; do awk 'BEGIN{OFS=" "} {match($1, /([0-9]+)-([0-9]+)/, a); diff = a[2] -a[1]; print $0, diff}' $i > ${i%.int.txt}".csv"; done
rm -r *.txt

for i in *.csv
do

tmp=$(mktemp)

awk 'BEGIN{OFS=" "} NR==1 && ($7 <= 2500) && ($8 >= ($13 - 2500))' $i > $tmp


if [[ -s $tmp ]]
then
cat $tmp >> insertii_conservate_rezultate.csv
else

echo $i >> fisiere_fara_rezultate.txt
fi
rm $tmp
done


mkdir insertii_specifice/
while read file; do   cp $file insertii_specifice/; done < fisiere_fara_rezultate.txt
rm fisiere_fara_rezultate.txt

cd insertii_specifice/
for file in *.csv; do

  awk '
  function abs(x) { return x < 0 ? -x : x }

  {
    lines[NR] = $0
    chrom[NR] = $2
    qstart[NR] = $7 + 0
    qend[NR] = $8 + 0
    sstart[NR] = $9 + 0
    send[NR] = $10 + 0
    qlen[NR] = $13 + 0
  }

  END {
 
    i = 1

    ins_len = qlen[i] - 6000  

    is_i_amonte = (qend[i] <= 3005)
    is_i_aval = (qstart[i] >= 3005)
    

    for (j = 2; j <= NR; j++) {

      if (chrom[i] != chrom[j]) continue

      is_j_amonte = (qend[j] <= 3005)
      is_j_aval = (qstart[j] >= 3005)

     
      if (is_i_amonte && is_j_aval &&
          qstart[j] - qend[i] >= ins_len - 20 &&
          qstar[j] - qend[i] <= ins_len + 20 ) {
         

        dist = abs(send[i] - sstart[j])
        if (dist <= 3000 - qend[i] + 50 ) {
          print lines[i] > "pereche_in_"FILENAME
          print lines[j] >> "pereche_in_"FILENAME
          exit
        }
      }

   
      if (is_i_aval && is_j_amonte &&
          qstart[i] - qend[j] >= ins_len - 20 &&
          qstart[i] - qend[j] <= ins_len + 20) {

        dist = abs(sstart[i] - send[j])
        if (dist <= 3000 - qend[j] + 50 ) {
          print lines[i] > "pereche_in_"FILENAME
          print lines[j] >> "pereche_in_"FILENAME
          exit
        }
      }
    }
  }
  ' "$file"
done

cat pereche* > insertii_specifice_rezultate.tsv
rm -r pereche*


while read name _ _ _ _ _ _ _ _ _ _ _ _; do id=${name%%:*}; echo $id; done < insertii_specifice_rezultate.tsv > lista_insertii.txt


for i in *.csv; do echo ${i%.csv}; done > lista_fisiere.txt


grep -Fvxf lista_insertii.txt lista_fisiere.txt > fara_informatii.txt

rm -r lista_insertii.txt lista_fisiere.txt
mkdir insertii_meh/
while read id; do cp ${id}*.csv insertii_meh/; done < fara_informatii.txt 

cd insertii_meh/

for file in *.csv; do
  awk '
  NR == 1 {
    best_line = $0;
    best_bitscore = $12;
    qstart = $7;
    qend = $8;

    if (qstart <= 500 && qend >= 3050) {
      keep_best = 1;
    } else {
      keep_best = 0;
    }
  }

  NR > 1 {
    diff = best_bitscore - $12;
    if (diff < 200) {
      keep_best = 0;
    }
  }

  END {
    if (keep_best) {
      print best_line > "fl_"FILENAME;
    }
  }
  ' "$file"
done

cat fl_* > insertii_conservate_un_fl_rezultate.csv
rm fl_* 
 

MAX_GAP=50
output="insertii_conservate_2hits_rezultate.csv"


> "$output"

for i in *.csv
do
  awk -v maxgap="$MAX_GAP" -v out="$output" '
  function abs(x) { return x < 0 ? -x : x }
  {
      qid = $1
      qstart = $7
      qend = $8
      sstart = $9
      send = $10
      qlen = $13

      alns[NR] = $0
      qs[NR] = qstart + 0
      qe[NR] = qend + 0
      ss[NR] = sstart + 0
      se[NR] = send + 0
      id[NR] = qid
      ql[NR] = qlen + 0
  }
  END {
      for (i = 1; i <= NR; i++) {
          for (j = 1; j <= NR; j++) {
              if (i == j) continue
              if (id[i] != id[j]) continue

              qlen = ql[i]
              min_q = (qs[i] < qs[j]) ? qs[i] : qs[j]
              max_q = (qe[i] > qe[j]) ? qe[i] : qe[j]
              gap_q = abs(qs[i] - qe[j])
              gap_s = abs(ss[i] - se[j])

              if (min_q <= 1 && max_q >= qlen && gap_q <= maxgap && gap_s <= maxgap) {
                  print alns[i] >> out
                  print alns[j] >> out
                  close(out)  # optional, pentru a închide fișierul după scriere
                  exit 0
              }
          }
      }
      
  }
  ' "$i"
done
cd ../../
mkdir results/
find . -type f -name "*_rezultate*" -exec cp -t results/ {} +
find . -type f -name "insertii" -exec cp -t results/ {} +
cd results/
for f in *.tsv; do
    mv "$f" "${f%.tsv}.csv"
done
awk -F ':' '{print $1}' *.csv | sort -u > identificate.txt

grep -Fvxf identificate.txt insertii > fara_informatii.txt


cp ../../all_ins.csv all_ins.csv


for f in insertii*.csv ; do mkdir ${f%.csv}; done


for f in *.csv; do dirname=${f%.csv}; if [ -d $dirname ]; then cp $f $dirname/; fi; done
##############################################################################################
cd insertii_conservate_rezultate/

awk '
FNR==NR {
  key = $1 "_" $2+1 "-" $3 ":" $6
  data[key] = $5 OFS $7
  next
}
{
  extra = ( $1 in data ) ? data[$1] : "NA\tNA"
  print $0, extra
}
' ../all_ins.csv insertii_conservate_rezultate.csv > insertii_conservate.csv


awk '{
    coord_ins = "NA"
    if ($13 < 6000 && $1 ~ /_2-/ && $9 > $10) {
        coord_ins = $9 - $14
    } else if ($13 < 6000 && $1 ~ /_2-/ && $9 < $10) {
        coord_ins = $9 + $14
    } else if ($13 > 6000 && $7 == 1 && $9 > $10) {
        coord_ins = $9 - 3000
    } else if ($13 > 6000 && $7 == 1 && $9 < $10) {
        coord_ins = $9 + 3000
    } else if ($13 > 6000 && $7 > 1 && $9 > $10) {
        coord_ins = $9 - (3000 - $7)
    } else if ($13 > 6000 && $7 > 1 && $9 < $10) {
        coord_ins = $9 + (3000 - $7)
    }
    print $0, coord_ins
}' insertii_conservate.csv > insertii_conservate_coord1.csv


awk '{

    coord_fin = "NA"

    # Calcul coord_fin
    delta = ($14 > $15) ? ($14 - $15) : ($15 - $14)
    if ($9 < $10) {
        coord_fin = $16 + delta
    } else if ($9 > $10) {
        coord_fin = $16 - delta
    }

    print $0, coord_fin
}' insertii_conservate_coord1.csv > insertii_conservate_coord.csv

rm insertii_conservate_coord1.csv

bash "$MATE_HOME/mate_associate.sh" insertii_conservate_coord.csv "$3" | sort -u > adnotari_conservate.txt
awk -F'\t' '{ids[$1]=1} END {for (id in ids) print id}' adnotari_conservate.txt > id_adnotari.txt
awk 'NR==FNR {id[$1]; next} !($1 in id)' id_adnotari.txt insertii_conservate_coord.csv > ins_conservate_neadnotate.txt
##############################################################################################
cd ../insertii_conservate_2hits_rezultate

awk '
FNR==NR {
  key = $1 "_" $2+1 "-" $3 ":" $6
  data[key] = $5 OFS $7
  next
}
{
  extra = ( $1 in data ) ? data[$1] : "NA\tNA"
  print $0, extra
}
' ../all_ins.csv insertii_conservate_2hits_rezultate.csv > insertii_conservate_2hits.csv


awk '($7 < 3000) && ($8 > (3000 + (($14 > $15) ? $14 - $15 : $15 - $14))) {
    offset = 3000 - $7;
    if ($9 < $10) {
        coord_ins = $9 + offset;
    } else {
        coord_ins = $9 - offset;
    }
    print $0, coord_ins;
}' insertii_conservate_2hits.csv > insertii_conservate_2hits_coord1.csv


awk '{

    coord_fin = "NA"

    # Calcul coord_fin
    delta = ($14 > $15) ? ($14 - $15) : ($15 - $14)
    if ($9 < $10) {
        coord_fin = $16 + delta
    } else if ($9 > $10) {
        coord_fin = $16 - delta
    }

    print $0, coord_fin
}' insertii_conservate_2hits_coord1.csv > insertii_conservate_2hits_coord.csv
rm insertii_conservate_2hits_coord1.csv


bash "$MATE_HOME/mate_associate.sh" insertii_conservate_2hits_coord.csv "$3" | sort -u > adnotari_conservate.txt
awk -F'\t' '{ids[$1]=1} END {for (id in ids) print id}' adnotari_2hits.txt > id_adnotari.txt
awk 'NR==FNR {id[$1]; next} !($1 in id)' id_adnotari.txt insertii_conservate_2hits_coord.csv > ins_conservate_2hits_neadnotate.txt

##############################################################################################
cd ../insertii_conservate_un_fl_rezultate

awk '
FNR==NR {
  key = $1 "_" $2+1 "-" $3 ":" $6
  data[key] = $5 OFS $7
  next
}
{
  extra = ( $1 in data ) ? data[$1] : "NA\tNA"
  print $0, extra
}
' ../all_ins.csv insertii_conservate_un_fl_rezultate.csv > insertii_conservate_fl.csv


awk '{
    coord_ins = "NA"
    if ($13 < 6000 && $1 ~ /_2-/ && $9 > $10) {
        coord_ins = $9 - $14
    } else if ($13 < 6000 && $1 ~ /_2-/ && $9 < $10) {
        coord_ins = $9 + $14
    } else if ($13 > 6000 && $7 == 1 && $9 > $10) {
        coord_ins = $9 - 3000
    } else if ($13 > 6000 && $7 == 1 && $9 < $10) {
        coord_ins = $9 + 3000
    } else if ($13 > 6000 && $7 > 1 && $9 > $10) {
        coord_ins = $9 - (3000 - $7)
    } else if ($13 > 6000 && $7 > 1 && $9 < $10) {
        coord_ins = $9 + (3000 - $7)
    }
    print $0, coord_ins
}' insertii_conservate_fl.csv > insertii_conservate_fl_coord1.csv


awk '{
coord_fin = "NA"
if ($9 < $10 ) {
coord_fin = $16 + ($8 - 3000)
} else if ($9 > $10) {
coord_fin = $16 - ($8 - 3000)
}
print $0, coord_fin
}' insertii_conservate_fl_coord1.csv > insertii_conservate_fl_coord.csv
rm insertii_conservate_fl_coord1.csv


bash "$MATE_HOME/mate_associate.sh" insertii_conservate_fl_coord.csv "$3" | sort -u > adnotari_conservate.txt
awk -F'\t' '{ids[$1]=1} END {for (id in ids) print id}' adnotari_fl.txt > id_adnotari.txt
awk 'NR==FNR {id[$1]; next} !($1 in id)' id_adnotari.txt insertii_conservate_fl_coord.csv > ins_conservate_fl_neadnotate.txt


cd ../insertii_specifice_rezultate

awk '{ 
  if ($8 <= 3005)
    print $1, $2, $10;
  else if ($7 > 3100)
    print $1, $2, $9;
}' insertii_specifice_rezultate.csv > coordonate_ins.csv


awk '{ 
  key = $1
  if (!(key in a)) {
    a[key] = $3   # salvează prima valoare a coloanei 3
    b[key] = $2   # salvează coloana 2 (e comună)
  } else {
    c[key] = $3   # salvează a doua valoare a coloanei 3
  }
}
END {
  for (k in a) {
    print k, b[k], a[k], c[k]
  }
}' coordonate_ins.csv > coord_insertii.csv
rm coordonate_ins.csv

cd $current_dir
done

mkdir results/
find . -type f -name coord_insertii.csv -exec cat {} + > results/specific_ins.csv
cd results/
awk '{
  split($1, id, "_")
  contig = id[1]  # doar partea de început ex: tig00000111
  chr = $2
  # sortăm coordonatele
  c1 = ($3 < $4) ? $3 : $4
  c2 = ($3 < $4) ? $4 : $3
  key = contig "_" chr "_" c1 "_" c2
  if (!(key in seen)) {
    seen[key] = 1
    print
  }
}' specific_ins.csv > specific_insertions.csv
rm specific_ins.csv
sed -i 's/ /,/g' specific_insertions.csv

cd ../
find . -type f -name adnotari_* -exec cat {} + > results/annotated_conserved_ins.csv
cd results/

awk '
{
    split($1, a, /_|[:-]/)
    len = a[3] - a[2]; if (len < 0) len = -len
    key = $2
    if (!(key in best) || len > best[key]) {
        best[key] = len
        line[key] = $0
    }
}
END {
    for (k in line) print line[k]
}' annotated_conserved_ins.csv > annotated_conserved_insertions.csv
rm annotated_conserved_ins.csv
sed -i 's/\t/,/g' annotated_conserved_insertions.csv

cd ../
find . -type f -name *_neadnotate.txt -exec cat {} + > results/unannotated_conserved_ins.csv

rm -r output* all_ins.csv *.bed *_blast.csv secv*

