#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/non_b_dna
encode_sign=all_chr_GQ
source ${public_script_dir}/public_common_encode.sh

