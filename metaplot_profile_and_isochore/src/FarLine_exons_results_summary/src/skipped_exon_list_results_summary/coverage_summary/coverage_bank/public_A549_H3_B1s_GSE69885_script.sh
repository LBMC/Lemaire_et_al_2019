#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_A549_H3_B1s_GSE69885
source ${public_script_dir}/public_common_simpler.sh

# maximum in heatmap
ymax=100

