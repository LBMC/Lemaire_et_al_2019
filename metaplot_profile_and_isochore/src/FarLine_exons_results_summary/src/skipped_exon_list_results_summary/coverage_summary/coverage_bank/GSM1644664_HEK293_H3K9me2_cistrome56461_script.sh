#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM1644664_HEK293_H3K9me2_cistrome56461
  source ${public_script_dir}/public_common_encode.sh
  
