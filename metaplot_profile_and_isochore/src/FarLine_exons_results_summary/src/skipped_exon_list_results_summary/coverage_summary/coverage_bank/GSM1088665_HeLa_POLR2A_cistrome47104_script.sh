#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/cistrome_PolII
  encode_sign=GSM1088665_HeLa_POLR2A_cistrome47104
  source ${public_script_dir}/public_common_encode.sh
  
