#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/histoneAc
encode_sign=GSM1059385_MCF-7_H3K122ac_33745
source ${public_script_dir}/public_common_encode.sh

