#!/bin/bash

exons_list2BED6_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${exons_list2BED6_dir}/exons_list2BED6.sh

exons_list2BED6_script () {
  local input_file=$1

  local out_dir=$( dirname ${input_file} )/bed6; mkdir -vp ${out_dir}
  local file_dir=${out_dir}/$( basename ${input_file} .tsv ).bed
  # echo "$file_dir" >&2
  exons_list2BED6 ${input_file} > ${file_dir}
}
