#!/bin/bash

exon_up_table=$1
exon_down_table=$2
FEATURE_TAB=$3
fig_dir=$4

exon_feature_plots_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for exon_table in ${exon_up_table} ${exon_down_table}; do
  Rscript ${exon_feature_plots_script_dir}/exon_feature_plots.r ${FEATURE_TAB} ${exon_table} ${fig_dir} --name $( basename ${exons_table} .tab )
done
