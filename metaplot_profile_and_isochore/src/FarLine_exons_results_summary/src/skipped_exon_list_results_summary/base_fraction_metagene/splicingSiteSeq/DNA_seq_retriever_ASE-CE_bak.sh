#!/bin/bash

DNA_seq_retriever_ASE_CE_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DNA_seq_retriever_ASE_CE_dir}/DNA_seq_retriever.sh

DNA_seq_retriever_ASE-CE () {
  ## process args
  local out_seq_dir=$1 #~ ./results/exons_sequences
  local out_ss_reg_dir=$2 #~ ./results/splicing_site_regions
  local metagene_param=$3
  local ASE_file=$4
  local CE_file=$5

  # source $CONF_FILE

  ## retrieve the DNA sequences
  # DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" "$window_size" "$in_exon" "$out_off" "$ASE_file" "$FASTA_DIR"
  # DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" "$window_size" "$in_exon" "$out_off" "$CE_file" "$FASTA_DIR"
  DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" $( echo ${metagene_param} | tr ',' ' ' ) "$ASE_file" "$FASTA_DIR"
  DNA_seq_retriever "$out_seq_dir" "$out_ss_reg_dir" $( echo ${metagene_param} | tr ',' ' ' ) "$CE_file" "$FASTA_DIR"
}
