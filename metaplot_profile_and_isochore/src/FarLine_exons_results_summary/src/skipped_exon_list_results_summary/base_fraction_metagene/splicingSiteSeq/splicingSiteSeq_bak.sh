#!/bin/bash

splicingSiteSeq_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${splicingSiteSeq_dir}/DNA_seq_retriever_uniq.sh

splicingSiteSeq () {
  ## Arguments
  local ssRegion_ext_list=$1
  local FASTA_DIR=$2
  local out_seq_dir=$3

  ## retrieve the corresponding DNA sequence
  DNA_seq_retriever_uniq $ssRegion_ext_list $FASTA_DIR ${out_seq_dir}/$( basename ${ssRegion_ext_list} .tsv ).seq
}
