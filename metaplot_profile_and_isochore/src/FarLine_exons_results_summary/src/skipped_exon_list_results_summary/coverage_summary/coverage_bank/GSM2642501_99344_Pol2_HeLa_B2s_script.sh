#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/ChIP_Pol2_15-11-18
  encode_sign=99344_Pol2_HeLa_B2s
  source ${public_script_dir}/public_common_encode.sh
  
