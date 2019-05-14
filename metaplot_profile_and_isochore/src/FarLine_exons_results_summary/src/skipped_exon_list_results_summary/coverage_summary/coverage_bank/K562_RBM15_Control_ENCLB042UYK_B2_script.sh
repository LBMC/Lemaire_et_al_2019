#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/RNA-Seq_Ctrl_SplicingLore
  encode_sign=K562_RBM15_Control_ENCLB042UYK_B2
  source ${public_script_dir}/public_common_encode.sh
  
