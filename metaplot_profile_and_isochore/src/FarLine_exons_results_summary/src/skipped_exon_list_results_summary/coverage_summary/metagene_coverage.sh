#!/bin/bash

metagene_coverage_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

metagene_coverage () {
  local bw_files=$1
  local cond=$2
  local rep=$3
  local bed6_files=$4
  local sign=$5
  local out_dir=$6
  local off_set=$7
  local comp_pair=$8
  local prefixes=$9

  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local hm_ymax=''
  if [[ "x${@}" == x*"--hm-ymax"* ]]; then
    hm_ymax="--hm-ymax $( echo ${@} | awk -F '--hm-ymax ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local mean_ymax=''
  if [[ "x${@}" == x*"--hm-ymax"* ]]; then
    mean_ymax="--ylims_mean $( echo ${@} | awk -F '--mean-ymax ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local median_ymax=''
  if [[ "x${@}" == x*"--hm-ymax"* ]]; then
    median_ymax="--ylims_median $( echo ${@} | awk -F '--median-ymax ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local sh_arg=''
  if [[ "x${@}" == x*"--scaleHarm"* ]]; then
    local sh_arg="--scaleHarm"
  fi

  local ref_chr_arg=''
  if [[ "x${@}" == x*'--ref-add-chr'* ]]; then
    local ref_chr_arg='--ref-add-chr'
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


  cmd="Rscript ${metagene_coverage_dir}/metagene_coverage.r \
  -bw ${bw_files} \
  -cond ${cond} \
  -rep ${rep} \
  -annot ${bed6_files} \
  -prefixes ${prefixes} \
  -sign ${sign} \
  -out_dir ${out_dir} \
  -off_set ${off_set} \
  -comp_pair ${comp_pair} \
  ${color_pallette} \
  ${hm_ymax} \
  ${mean_ymax} \
  ${median_ymax} \
  ${sh_arg} \
  ${ref_chr_arg} \
  ${various_length} \
  ${side_cut} \
  ${meta_2pos_ext} \
  ${gene_concat_arg} \
  ;"
  echo "$cmd"
  #~ eval "$cmd"

  qsub_script=${out_dir}/qsub_script.sh
  less ${metagene_coverage_dir}/run_metagene_coverage.sh | sed s%'ppp_cmd_ppp'%"$cmd"%  | sed s%'ppp_logs_[erout]*_dir_ppp'%"${out_dir}"% > $qsub_script
  cmd="qsub $qsub_script"
  echo "$cmd"
  eval "$cmd"

}
