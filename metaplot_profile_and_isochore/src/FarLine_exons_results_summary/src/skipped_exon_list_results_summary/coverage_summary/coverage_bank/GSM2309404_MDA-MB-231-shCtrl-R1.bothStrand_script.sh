#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/run_on
  encode_sign=GSM2309404_MDA-MB-231-shCtrl-R1.bothStrand
  source ${public_script_dir}/public_common_encode.sh
  
