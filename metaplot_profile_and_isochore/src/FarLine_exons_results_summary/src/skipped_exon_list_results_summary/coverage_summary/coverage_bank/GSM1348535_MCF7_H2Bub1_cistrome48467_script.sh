#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM1348535_MCF7_H2Bub1_cistrome48467
  source ${public_script_dir}/public_common_encode.sh
  
