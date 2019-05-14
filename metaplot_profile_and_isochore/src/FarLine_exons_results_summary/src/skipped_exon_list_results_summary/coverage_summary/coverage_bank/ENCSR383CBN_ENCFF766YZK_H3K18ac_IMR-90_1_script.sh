#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K18ac
encode_sign=ENCSR383CBN_ENCFF766YZK_H3K18ac_IMR-90_1
source ${public_script_dir}/public_common_encode.sh

