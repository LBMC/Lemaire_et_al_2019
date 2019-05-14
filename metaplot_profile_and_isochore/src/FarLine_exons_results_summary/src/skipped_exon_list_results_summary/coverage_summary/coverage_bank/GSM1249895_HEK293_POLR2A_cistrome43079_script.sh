#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/cistrome_PolII
  encode_sign=GSM1249895_HEK293_POLR2A_cistrome43079
  source ${public_script_dir}/public_common_encode.sh
  
