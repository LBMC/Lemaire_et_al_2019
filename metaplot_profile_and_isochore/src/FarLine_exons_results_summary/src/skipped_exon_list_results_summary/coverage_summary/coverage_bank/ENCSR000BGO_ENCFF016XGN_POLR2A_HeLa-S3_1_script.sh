#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/POLR2A
encode_sign=ENCSR000BGO_ENCFF016XGN_POLR2A_HeLa-S3_1
source ${public_script_dir}/public_common_encode.sh

