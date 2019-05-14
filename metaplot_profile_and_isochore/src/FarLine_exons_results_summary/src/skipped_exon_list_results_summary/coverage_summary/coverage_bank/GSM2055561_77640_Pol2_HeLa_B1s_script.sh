#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/ChIP_Pol2_15-11-18
  encode_sign=77640_Pol2_HeLa_B1s
  source ${public_script_dir}/public_common_encode.sh
  
