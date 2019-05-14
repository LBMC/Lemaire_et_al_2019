#!/bin/bash

for xxx in $( ls SFs/ ); do
  if [[ "x$( head -1 SFs/$xxx)" != x'exon_gene_pos' ]]; then
    echo "exon_gene_pos" > SFs_seb/$xxx
    less SFs/${xxx} >> SFs_seb/$xxx
  else
    less SFs/${xxx} > SFs_seb/$xxx
  fi
done

bash nb_exon_list.sh

bash GCbias_exonList2FLtab.sh
bash GCbias_nbExonList2FLtab.sh

####
