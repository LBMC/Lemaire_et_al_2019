#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/RNA-Seq_Ctrl_SplicingLore
  encode_sign=293T_ESRP2_neg_ctrl
  source ${public_script_dir}/public_common_encode.sh
  
