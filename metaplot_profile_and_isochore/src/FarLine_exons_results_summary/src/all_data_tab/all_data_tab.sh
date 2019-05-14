#!/bin/bash

out_dir=$1
exons_list_name=$2
all_data_tab=$3
GENE_FDB_TAB=$4
EXON_FDB_TAB=$5

all_data_tab_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


deltaPsi_thresh='0.10'
pval_thresh='0.05'

exons_up_table=${out_dir}/exon_up_${exons_list_name}.tab
bash ${all_data_tab_dir}/table_exons_builder.sh ${all_data_tab} 'up' ${GENE_FDB_TAB} ${EXON_FDB_TAB} ${deltaPsi_thresh} ${pval_thresh} > ${exons_up_table}

exons_down_table=${out_dir}/exon_down_${exons_list_name}.tab
bash ${all_data_tab_dir}/table_exons_builder.sh ${all_data_tab} 'down' ${GENE_FDB_TAB} ${EXON_FDB_TAB} ${deltaPsi_thresh} ${pval_thresh} > ${exons_down_table}
