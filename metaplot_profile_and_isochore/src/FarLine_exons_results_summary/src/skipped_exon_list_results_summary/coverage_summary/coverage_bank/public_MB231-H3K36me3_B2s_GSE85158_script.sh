#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/bw_files/public_MB231-H3K36me3_B2s_GSE85158
source ${public_script_dir}/public_common.sh

# maximum in heatmap
ymax=38.712288183956545

