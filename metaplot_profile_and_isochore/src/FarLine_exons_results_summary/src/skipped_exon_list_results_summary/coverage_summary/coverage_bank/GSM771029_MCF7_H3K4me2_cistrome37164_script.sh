#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM771029_MCF7_H3K4me2_cistrome37164
  source ${public_script_dir}/public_common_encode.sh
  
