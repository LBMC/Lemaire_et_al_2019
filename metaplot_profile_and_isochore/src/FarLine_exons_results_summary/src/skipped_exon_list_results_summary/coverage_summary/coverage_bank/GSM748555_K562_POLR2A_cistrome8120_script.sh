#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/cistrome_PolII
  encode_sign=GSM748555_K562_POLR2A_cistrome8120
  source ${public_script_dir}/public_common_encode.sh
  
