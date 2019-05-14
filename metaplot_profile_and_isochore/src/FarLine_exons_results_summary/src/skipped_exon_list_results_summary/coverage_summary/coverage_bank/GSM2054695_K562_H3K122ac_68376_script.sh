#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/histoneAc
encode_sign=GSM2054695_K562_H3K122ac_68376
source ${public_script_dir}/public_common_encode.sh

