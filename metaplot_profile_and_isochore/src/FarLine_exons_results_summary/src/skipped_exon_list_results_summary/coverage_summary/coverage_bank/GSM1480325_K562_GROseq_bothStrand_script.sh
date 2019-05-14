#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/run_on
  encode_sign=GSM1480325_K562_GROseq_bothStrand
  source ${public_script_dir}/public_common_encode.sh
  
