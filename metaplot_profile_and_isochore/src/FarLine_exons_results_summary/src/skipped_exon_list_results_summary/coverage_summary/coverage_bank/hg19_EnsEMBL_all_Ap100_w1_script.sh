#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/GCp100
encode_sign=hg19_EnsEMBL_all_Ap100_w1
source ${public_script_dir}/public_common_encode.sh

