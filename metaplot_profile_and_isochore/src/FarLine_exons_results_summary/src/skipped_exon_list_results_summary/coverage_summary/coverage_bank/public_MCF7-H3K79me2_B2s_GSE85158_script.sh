#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_MCF7-H3K79me2_B2s_GSE85158
source ${public_script_dir}/public_common.sh

# maximum in heatmap
ymax=89.99821963733847

