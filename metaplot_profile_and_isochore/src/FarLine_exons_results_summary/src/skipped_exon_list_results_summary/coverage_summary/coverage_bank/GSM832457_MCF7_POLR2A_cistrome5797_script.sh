#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/cistrome_PolII
  encode_sign=GSM832457_MCF7_POLR2A_cistrome5797
  source ${public_script_dir}/public_common_encode.sh
  
