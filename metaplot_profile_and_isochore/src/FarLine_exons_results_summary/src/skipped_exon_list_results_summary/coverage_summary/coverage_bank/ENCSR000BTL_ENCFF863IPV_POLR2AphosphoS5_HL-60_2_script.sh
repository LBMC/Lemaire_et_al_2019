#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/POLR2AphosphoS5
encode_sign=ENCSR000BTL_ENCFF863IPV_POLR2AphosphoS5_HL-60_2
source ${public_script_dir}/public_common_encode.sh

