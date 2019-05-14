#!/bin/bash

kmeans_clustering_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

kmeans_clustering () {
  local bw_files=$1
  local cond=$2
  local rep=$3
  local bed6_file=$4
  local sign=$5
  local out_dir=$6
  local off_set=$7
  local comp_pair=$8
  local nbr_cluster=$9
  local min_max_signal=${10}

  local hm_ymax=''
  if [[ "x${@}" == x*"--hm-ymax"* ]]; then
    hm_ymax="--hm-ymax $( echo ${@} | awk -F '--hm-ymax ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi


  Rscript ${kmeans_clustering_dir}/kmeans_clustering.r \
  -bw ${bw_files} \
  -cond ${cond} \
  -rep ${rep} \
  -annot ${bed6_file} \
  -sign ${sign} \
  -out_dir ${out_dir} \
  -off_set ${off_set} \
  -comp_pair ${comp_pair} \
  -nbr_cluster ${nbr_cluster} \
  -min_max_signal ${min_max_signal} \
  ${hm_ymax} \
  ;
  # -scaleHarm \
}
