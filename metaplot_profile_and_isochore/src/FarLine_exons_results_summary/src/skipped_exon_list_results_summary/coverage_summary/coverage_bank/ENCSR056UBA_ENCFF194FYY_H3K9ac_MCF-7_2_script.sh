#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/encode_histoneMarks/H3K9ac
encode_sign=ENCSR056UBA_ENCFF194FYY_H3K9ac_MCF-7_2
source ${public_script_dir}/public_common_encode.sh

