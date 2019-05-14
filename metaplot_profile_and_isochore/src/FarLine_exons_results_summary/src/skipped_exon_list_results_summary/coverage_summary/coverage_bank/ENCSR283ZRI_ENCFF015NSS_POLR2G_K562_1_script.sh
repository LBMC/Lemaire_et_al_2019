#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/POLR2G
encode_sign=ENCSR283ZRI_ENCFF015NSS_POLR2G_K562_1
source ${public_script_dir}/public_common_encode.sh

