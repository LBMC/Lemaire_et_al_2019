#!/bin/bash

methylation_metagene_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${methylation_metagene_script_dir}/methylation_metagene.sh

methylation_metagene_script () {
  local seq_files=$1
  local meth_files=$2
  local fig_dir=$3
  local SS=$4
  local metagene_param=$5
  local names=$6
  local expSeq_cond=$7
  local expSeq_rep=$8
  local comp_pair=$9
  local min_reads=${10}

  local NBRCORES=1
  if [[ "x$@" == x*"--nbr-cores"* ]]; then
    NBRCORES=$( echo "$@" | awk -F '--nbr-cores' '{ print $NF }' | cut -d ' ' -f 2 )
  fi

  local fig_prefix='test'
  if [[ "x$@" == x*"--fig-prefix"* ]]; then
    local fig_prefix=$( echo "$@" | awk -F '--fig-prefix' '{ print $NF }' | cut -d ' ' -f 2 )
  fi

  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  cmd="
  methylation_metagene \
    ${seq_files} \
    ${meth_files} \
    ${names} \
    ${SS} \
    $fig_dir \
    $( echo ${metagene_param} | tr ',' ' ' ) \
    ${expSeq_cond} \
    ${expSeq_rep} \
    ${comp_pair} \
    ${min_reads} \
    --fig-prefix ${fig_prefix} \
    --nbr-cores $NBRCORES \
    ${color_pallette}
    "
    # echo "$cmd"
    eval "$cmd"
}
