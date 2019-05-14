#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/cistrome_histoneUb
encode_sign=GSE94973_U2OS_H2AK119ub_1
source ${public_script_dir}/public_common_encode.sh

