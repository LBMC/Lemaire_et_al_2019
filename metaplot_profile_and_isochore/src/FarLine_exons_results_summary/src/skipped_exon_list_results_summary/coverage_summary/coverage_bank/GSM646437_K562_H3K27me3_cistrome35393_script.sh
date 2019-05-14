#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM646437_K562_H3K27me3_cistrome35393
  source ${public_script_dir}/public_common_encode.sh
  
