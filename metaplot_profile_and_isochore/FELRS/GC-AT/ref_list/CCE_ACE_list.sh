#!/bin/bash

tab_db=~/analyses_SEBASTIEN/data/tab_files/Nicolas_db/hsapiens_exonsstatus_improved_full.csv

exon_type_col="$( head -1 $tab_db | tr '\t' '\n' | grep -En '^exon_types$' | cut -d ':' -f 1 )"

for extype in CCE ACE; do
  out_file=./${extype}.txt
  echo "exon_gene_pos" > ${out_file}
  tail -n +2 $tab_db | awk -F '\t' -v extype=$extype -v etc=$exon_type_col '{ if ( $etc == extype ) print $1"_"$2 }' >> ${out_file}
done


####
