#!/bin/bash

pptx_builder_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

fig_dir=$1
pptx_name=$2
exon_list_name=$3
res_dir=$4
CONF_FILE=$5

python3 ${pptx_builder_dir}/pptx_builder.py ${fig_dir} ${pptx_name} ${exon_list_name} ${res_dir}/exons_lists ${CONF_FILE}
