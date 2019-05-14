#!/bin/bash

splicingSiteSeq_ASE_CE_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${splicingSiteSeq_ASE_CE_dir}/DNA_seq_retriever_uniq.sh

splicingSiteSeq_ASE-CE () {
  ## Arguments
  local ssRegion_ext_dir=$1
  local FASTA_DIR=$2 #~ ~/analyses_SEBASTIEN/data/genome/hg19/EnsEMBL_all
  local out_seq_dir=$3 #~ ./results/exons_sequences
  local SS=$4

  # rebuild path to ASE and CE splicing site region lists
  local ASE_list=${ssRegion_ext_dir}/ASE_FDB_${SS}_ext.tsv
  local CE_list=${ssRegion_ext_dir}/CE_FDB_${SS}_ext.tsv

  ## recover the sequences in a region around each splicing site, with extension to compute percentage in a window centered on each base
  # DNA_seq_retriever_ASE-CE "$out_seq_dir" "$out_ss_reg_dir" "$metagene_param" "$ASE_list" "$CE_list" "$FASTA_DIR"
  local ssReg_list=''
  for ssReg_list in ${ASE_list} ${CE_list}; do
    DNA_seq_retriever_uniq $ssReg_list $FASTA_DIR ${out_seq_dir}/$( basename ${ssReg_list} .tsv ).seq
  done
}
