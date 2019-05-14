#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_AML-Gfi-36S_H3_B2s_GSE71254
source ${public_script_dir}/public_common_simpler.sh

# maximum in heatmap
ymax=100

