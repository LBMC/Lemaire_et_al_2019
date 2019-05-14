#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_MCF7-H3K4me1_B1s_GSE85158
source ${public_script_dir}/public_common.sh

# maximum in heatmap
ymax=40.134076568526986

