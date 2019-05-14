#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files
  encode_sign=GSM2135493_GSE80728_MCF7_Pol2_dl
  source ${public_script_dir}/public_common_encode.sh
  
