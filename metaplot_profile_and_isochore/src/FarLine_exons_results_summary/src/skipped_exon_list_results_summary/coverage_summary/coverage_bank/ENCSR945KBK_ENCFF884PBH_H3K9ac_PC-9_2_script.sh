#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K9ac
encode_sign=ENCSR945KBK_ENCFF884PBH_H3K9ac_PC-9_2
source ${public_script_dir}/public_common_encode.sh

