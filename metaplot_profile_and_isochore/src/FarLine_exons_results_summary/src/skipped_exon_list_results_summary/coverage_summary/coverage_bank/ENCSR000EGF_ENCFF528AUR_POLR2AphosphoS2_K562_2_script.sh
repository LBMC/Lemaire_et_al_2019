#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/POLR2AphosphoS2
encode_sign=ENCSR000EGF_ENCFF528AUR_POLR2AphosphoS2_K562_2
source ${public_script_dir}/public_common_encode.sh

