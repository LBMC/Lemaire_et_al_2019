#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/Nuc_20-11-18
  encode_sign=GSM971951_GSE39579_MNase_HEK293_B1s
  source ${public_script_dir}/public_common_encode.sh
  