#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K27ac
encode_sign=ENCSR007YOT_ENCFF576XLQ_H3K27ac_GM23248_2
source ${public_script_dir}/public_common_encode.sh

