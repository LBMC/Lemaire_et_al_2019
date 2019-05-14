#!/bin/bash

seqFile2Fasta () {
  local seq_file_list=$1
  local out_fasta_list=$2

  local seq_file_array_spec=( $( echo "${seq_file_list}" | tr ',' ' ' ) )
  local out_fasta_array_spec=( $( echo "${out_fasta_list}" | tr ',' ' ' ) )

  for (( iii=0; iii<${#seq_file_array_spec[*]}; iii++ )); do
    cmd="
    tail -n +2 ${seq_file_array_spec[$iii]} | awk '{ print \">\"\$1\":\"\$2\"\n\"\$3 }' > ${out_fasta_array_spec[$iii]}
    "

    # echo "$cmd"
    eval "$cmd"
  done

}
