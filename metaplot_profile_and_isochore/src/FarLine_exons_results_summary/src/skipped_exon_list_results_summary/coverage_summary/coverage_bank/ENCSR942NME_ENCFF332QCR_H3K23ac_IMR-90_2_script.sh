#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K23ac
encode_sign=ENCSR942NME_ENCFF332QCR_H3K23ac_IMR-90_2
source ${public_script_dir}/public_common_encode.sh

