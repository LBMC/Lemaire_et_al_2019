#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_HeLa_H3_B1s_GSE60526
source ${public_script_dir}/public_common_simpler.sh

# maximum in heatmap
ymax=100
