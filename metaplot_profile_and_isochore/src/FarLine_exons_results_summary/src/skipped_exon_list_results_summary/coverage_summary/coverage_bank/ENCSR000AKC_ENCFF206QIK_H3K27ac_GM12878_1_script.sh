#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K27ac
encode_sign=ENCSR000AKC_ENCFF206QIK_H3K27ac_GM12878_1
source ${public_script_dir}/public_common_encode.sh

