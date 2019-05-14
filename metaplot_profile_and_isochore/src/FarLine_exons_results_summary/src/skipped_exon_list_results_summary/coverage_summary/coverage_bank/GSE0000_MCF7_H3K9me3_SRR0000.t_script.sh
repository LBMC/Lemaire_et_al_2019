#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  expSeq_bw_dir=${data_dir}/bw_files/ChromLore
  encode_sign=GSE0000_MCF7_H3K9me3_SRR0000.t
  source ${public_script_dir}/public_common_encode.sh
  
