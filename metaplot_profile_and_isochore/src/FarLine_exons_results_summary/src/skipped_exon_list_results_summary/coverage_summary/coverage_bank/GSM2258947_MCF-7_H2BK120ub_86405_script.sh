#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/histoneUb
encode_sign=GSM2258947_MCF-7_H2BK120ub_86405
source ${public_script_dir}/public_common_encode.sh

