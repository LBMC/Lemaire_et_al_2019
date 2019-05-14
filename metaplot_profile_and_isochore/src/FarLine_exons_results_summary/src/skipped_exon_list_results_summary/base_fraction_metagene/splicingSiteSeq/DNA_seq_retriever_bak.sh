#!/bin/bash

DNA_seq_retriever_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DNA_seq_retriever_dir}/DNA_seq_retriever_uniq.sh

DNA_seq_retriever () {
  ## process args
  local out_seq_dir=$1 #~ ./results/exons_sequences
  local window_size=$3 #~ 20
  local in_exon=$4 #~ 50
  local in_off=$(($in_exon-1))
  local out_off=$5 #~ 100
  local file_exons_list=$6 #~ ${ex_list_dir}/exon_list.tsv
  local FASTA_DIR=$7

  ## retrive the sequences for list of exons, list of ASE and list of CE
  DNA_seq_retriever_uniq $reg_3SS_ext_file $FASTA_DIR $out_3SS_fa
  DNA_seq_retriever_uniq $reg_5SS_ext_file $FASTA_DIR $out_5SS_fa

}
