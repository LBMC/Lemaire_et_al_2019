#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM335958_HeLa_H3.3_H2A.Z_cistrome1234
  source ${public_script_dir}/public_common_encode.sh
  
