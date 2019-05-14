#!/bin/bash

metagene_diffCoverage_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

metagene_diffCoverage () {
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

  local sh_arg=''
  if [[ "x${@}" == x*"--scaleHarm"* ]]; then
    sh_arg="--scaleHarm"
  fi

  Rscript ${metagene_diffCoverage_dir}/metagene_diffCoverage.r \
  -bw ${bw_files} \
  -cond ${cond} \
  -rep ${rep} \
  -annot ${bed6_files} \
  -prefixes ${prefixes} \
  -sign ${sign} \
  -out_dir ${out_dir} \
  -off_set ${off_set} \
  -comp_pair ${comp_pair} \
  -ylims_mean NULL,NULL \
  -ylims_median NULL,NULL \
  ${color_pallette} \
  ${sh_arg} \
  ;
  # -scaleHarm \
}
