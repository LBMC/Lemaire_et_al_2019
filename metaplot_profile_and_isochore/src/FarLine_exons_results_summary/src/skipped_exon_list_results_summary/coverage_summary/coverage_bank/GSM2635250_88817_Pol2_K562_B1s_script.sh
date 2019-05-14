#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/ChIP_Pol2_15-11-18
  encode_sign=88817_Pol2_K562_B1s
  source ${public_script_dir}/public_common_encode.sh
  
