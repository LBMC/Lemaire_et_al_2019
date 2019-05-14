#!/bin/bash

metagene_plots_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${metagene_plots_script_dir}/metagene_plots.sh

metagene_plots_script () {
  local seq_files=$1
  local fig_dir=$2
  local SS=$3
  local metagene_param=$4
  local names=$5

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
  metagene_plots \
    ${seq_files} \
    ${names} \
    ${SS} \
    $fig_dir \
    $( echo ${metagene_param} | tr ',' ' ' ) \
    --fig-prefix ${fig_prefix} \
    --nbr-cores $NBRCORES \
    ${color_pallette}
    "
    # echo "$cmd"
    eval "$cmd"
}
