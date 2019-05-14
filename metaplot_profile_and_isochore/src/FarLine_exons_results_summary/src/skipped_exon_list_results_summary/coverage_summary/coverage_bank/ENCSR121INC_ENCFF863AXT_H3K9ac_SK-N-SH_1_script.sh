#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K9ac
encode_sign=ENCSR121INC_ENCFF863AXT_H3K9ac_SK-N-SH_1
source ${public_script_dir}/public_common_encode.sh

