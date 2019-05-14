#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/cistrome_histoneUb
encode_sign=GSE85158_MB436-H2BK120ub
source ${public_script_dir}/public_common_encode.sh

