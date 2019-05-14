#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/BIGWIG_hg38
  encode_sign=GSM916107_MCF7_H3K36me3_cistrome38326
  source ${public_script_dir}/public_common_encode.sh
  
