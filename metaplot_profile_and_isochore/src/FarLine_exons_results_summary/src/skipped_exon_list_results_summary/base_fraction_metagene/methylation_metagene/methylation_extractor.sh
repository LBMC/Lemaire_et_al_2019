#!/bin/bash

methylation_extractor_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

methylation_extractor () {
  local bisulfite_bam=$1
  local seq_file_list=$2
  local meth_file_list=$3

  local seq_files_all=( $( echo ${seq_file_list} | tr ',' ' ' ) )
  local meth_file_array=( $( echo ${meth_file_list} | tr ',' ' ' ) )

  for (( idx=0; idx<${#meth_file_array[*]}; idx++ )); do
    local out_meth_file=${meth_file_array[$idx]}
    local seq_file=${seq_files_all[$idx]}

    cmd="
    python3 ${methylation_extractor_dir}/methylation_extractor.py ${bisulfite_bam} ${seq_file} > ${out_meth_file}
    "
    echo "$cmd"
    eval "$cmd"
  done

}
