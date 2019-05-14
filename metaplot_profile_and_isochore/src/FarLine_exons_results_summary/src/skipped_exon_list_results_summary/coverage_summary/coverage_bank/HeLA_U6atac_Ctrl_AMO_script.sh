#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/RNA-Seq_Ctrl_SplicingLore
  encode_sign=HeLA_U6atac_Ctrl_AMO
  source ${public_script_dir}/public_common_encode.sh
  
