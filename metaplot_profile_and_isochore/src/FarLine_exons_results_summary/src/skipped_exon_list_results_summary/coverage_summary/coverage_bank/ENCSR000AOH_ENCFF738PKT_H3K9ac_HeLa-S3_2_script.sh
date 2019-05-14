#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K9ac
encode_sign=ENCSR000AOH_ENCFF738PKT_H3K9ac_HeLa-S3_2
source ${public_script_dir}/public_common_encode.sh

