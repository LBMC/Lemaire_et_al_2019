#!/bin/bash

all_data_tab=$1
exon_type=$2
GENE_FDB_TAB=$3
EXON_FDB_TAB=$4
abs_deltaPsi_thresh=$5
pval_thresh=$6

table_exons_builder_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## Some functions
source ${table_exons_builder_dir}/col_nbr.sh

filtre () {
  exon_type=$1
  all_txt="$( more "${2:-/dev/stdin}" )"
  deltaPSI_col_nbr="$( echo "${all_txt}" | col_nbr 'deltaPSI' )"

  echo "${all_txt}" | head -1

  if [[ "x$exon_type" == x"up" ]]; then
    echo "${all_txt}" | tail -n +2 | awk -v delta_col=${deltaPSI_col_nbr} '( $delta_col > 0 )'
  elif [[ "x$exon_type" == x"down" ]]; then
    echo "${all_txt}" | tail -n +2 | awk -v delta_col=${deltaPSI_col_nbr} '( $delta_col < 0 )'
  elif [[ "x$exon_type" == x"all" ]]; then
    echo "${all_txt}" | tail -n +2 
  fi
}


## determine the column index of for id_gene
id_gene_col_nbr="$( col_nbr 'id_gene' ${all_data_tab}  )" #~15"
id_exon_col_nbr="$( col_nbr 'exon_skipped' ${all_data_tab}  )" #~15"

## write the headers
echo -e "id_event\tid_exon\texons_flanquants\tcoordonnees\tid_gene\tgene_symbol\texon_pos"

## add the FasterDB exon id to the 'all_data_tab' table, only for selected exon based on statistical test and deltaPsi

join -j 1 \
  <(
    Rscript ${table_exons_builder_dir}/significant_exon_selection.r ${all_data_tab} --pval ${pval_thresh} --delta ${abs_deltaPsi_thresh} | filtre ${exon_type} | tail -n +2 | awk -v exon_col=$id_exon_col_nbr -v gene_col=$id_gene_col_nbr '{ OFS="\t"; print $gene_col"_"$exon_col,$1,$5,$4,$gene_col,$2,$exon_col }' | sort -k1,1
  ) \
  <(
    tail -n +2 $EXON_FDB_TAB | awk '{ OFS="\t"; print $2"_"$3,$1 }' | sort -k1,1
  ) | awk '{ OFS="\t"; print $2,$8,$3,$4,$5,$6,$7 }'
