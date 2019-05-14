#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/ChIP_Pol2_15-11-18
  encode_sign=75792_Pol2S2P_K562_B2s
  source ${public_script_dir}/public_common_encode.sh
  
