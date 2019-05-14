#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/POLR2B
encode_sign=ENCSR325RLL_ENCFF545XLF_POLR2B_K562_2
source ${public_script_dir}/public_common_encode.sh

