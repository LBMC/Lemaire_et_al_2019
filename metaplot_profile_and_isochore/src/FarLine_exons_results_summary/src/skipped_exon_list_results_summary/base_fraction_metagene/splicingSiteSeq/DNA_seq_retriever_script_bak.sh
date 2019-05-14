#!/bin/bash

DNA_seq_retriever_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DNA_seq_retriever_script_dir}/DNA_seq_retriever.sh

## process args
out_seq_dir=$1 #~ ./results/exons_sequences
out_ss_reg_dir=$2 #~ ./results/splicing_site_regions
CONF_FILE=$3
exons_up_list=$4
exons_down_list=$5

source $CONF_FILE

DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" "$window_size" "$in_exon" "$out_off" "$exons_up_list" "$FASTA_DIR"
DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" "$window_size" "$in_exon" "$out_off" "$exons_down_list" "$FASTA_DIR"
