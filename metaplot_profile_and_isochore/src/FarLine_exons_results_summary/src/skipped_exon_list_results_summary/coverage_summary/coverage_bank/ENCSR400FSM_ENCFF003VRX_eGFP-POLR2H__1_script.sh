#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
expSeq_bw_dir=${data_dir}/bw_files/encode_Pol2/eGFP-POLR2H
encode_sign=ENCSR400FSM_ENCFF003VRX_eGFP-POLR2H__1
source ${public_script_dir}/public_common_encode.sh

