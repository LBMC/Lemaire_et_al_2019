#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/histoneUb
encode_sign=GSM818830_HeLa_H2BK120ub_5787
source ${public_script_dir}/public_common_encode.sh

