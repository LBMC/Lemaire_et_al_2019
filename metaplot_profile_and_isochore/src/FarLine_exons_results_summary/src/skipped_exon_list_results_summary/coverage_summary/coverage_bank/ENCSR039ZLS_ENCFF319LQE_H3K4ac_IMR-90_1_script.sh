#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K4ac
encode_sign=ENCSR039ZLS_ENCFF319LQE_H3K4ac_IMR-90_1
source ${public_script_dir}/public_common_encode.sh

