#!/bin/bash

DNA_seq_retriever_uniq_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

DNA_seq_retriever_uniq () {
  local reg_coord_file_list=$1
  local FASTA_DIR=$2
  local out_file_list=$3

  python3 ${DNA_seq_retriever_uniq_dir}/list_coord_seq_retriever.py ${FASTA_DIR} $reg_coord_file_list $out_file_list

  # remove the duplicates
  for out_file in $( echo $out_file_list | tr ',' ' ' ); do
    head -1 $out_file > ${out_file}_tmp
    tail -n +2 $out_file  | sort | uniq >> ${out_file}_tmp
    mv ${out_file}_tmp $out_file
  done
}
