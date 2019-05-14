#!/bin/bash

coverage_means_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

coverage_means () {
  ## process arguments
  local expSeq_bw_files=$1
  local bed6_files=$2 
  local prefixes=$3
  local out_dir=$4
  local sign=$5

  mkdir -vp $out_dir

  flag='--color-pallette'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local sum_cov_chr_arg=''
  flag='--ref-add-chr'
  if [[ "x${@}" == x*"${flag}"* ]]; then
      local sum_cov_chr_arg='--add-chr 1'
  fi

  local various_length=''
  flag='--various-length'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local various_length="--various-length $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"

    local side_cut=''
    flag='--side-cut'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local side_cut="--side-cut $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi

    local meta_2pos_ext=''
    flag='--ext-up'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local meta_2pos_ext="${meta_2pos_ext} --ext-up $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi
    flag='--ext-dw'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local meta_2pos_ext="${meta_2pos_ext} --ext-dw $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi
  fi
  # local meta_2pos_ext='--ext-up 200 --ext-dw 200'

  local gene_concat_arg=''
  if [[ "x${@}" == x*"--gene-concat"* ]]; then
    local gene_concat_arg="--gene-concat"
  fi


  ebf_array=( $( echo ${expSeq_bw_files} | tr ',' ' ' ) )

  bed6_array=( $( echo ${bed6_files} | tr ',' ' ' ) )
  prefixes_array=( $( echo ${prefixes} | tr ',' ' ' ) )

  for cov_idx in $( seq 1 ${#ebf_array[*]} ); do
    real_cov_idx=$((cov_idx-1))
    mean_files_array=()
    for bed_idx in $( seq 1 ${#bed6_array[*]} ); do
      real_bed_idx=$((bed_idx-1))
      out_file=${out_dir}/${prefixes_array[$real_bed_idx]}_${sign}_means.tab

      mean_files_array=( ${mean_files_array[*]} ${out_file} )
      #~ cmd="python3 ${coverage_means_dir}/metagene_coverage_functions/retrieve_bw_values_mean.py ${ebf_array[$real_cov_idx]} ${bed6_array[$real_bed_idx]} ${out_file} ${sum_cov_chr_arg} --with-ids"
      cmd="Rscript ${coverage_means_dir}/coverage_means.r ${ebf_array[$real_cov_idx]} ${bed6_array[$real_bed_idx]} ${out_file} ${sum_cov_chr_arg} ${various_length} ${side_cut} ${meta_2pos_ext} ${gene_concat_arg}" # --with-ids
      echo "$cmd"
      eval "$cmd" #~ >> ${out_file}
    done

    for cov_idx in $( seq 1 ${#ebf_array[*]} ); do
      mean_files="$( echo ${mean_files_array[*]} | tr ' ' ',' )"
      cmd="Rscript ${coverage_means_dir}/coverage_means_bxpl_stat.r ${mean_files} ${prefixes} ${out_dir} ${sign} ${color_pallette}"
      echo "$cmd"
      eval "$cmd"
    done

  done
}
