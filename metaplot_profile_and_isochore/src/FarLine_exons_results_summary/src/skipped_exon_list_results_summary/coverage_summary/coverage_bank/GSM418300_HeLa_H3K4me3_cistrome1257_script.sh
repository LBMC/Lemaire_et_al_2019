#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM418300_HeLa_H3K4me3_cistrome1257
  source ${public_script_dir}/public_common_encode.sh
  
