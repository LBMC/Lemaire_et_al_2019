#!/bin/bash

coverage_summary_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${coverage_summary_script_dir}/coverage_summary.sh

coverage_summary_script () {
  local input_prefixes=$1
  local bed_dir=$2
  local fig_dir=$3
  local CONF_FILE=$4
  local data_dir=$5
  # local prefixes_ref=$6
  # local cov_in_exon=$4
  # local cov_out_off=$5
  # local data_dir=$6

  local prefixes_ref=''
  flag='--ex-tab-ref'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local prefixes_ref="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    local prefixes_ref_arg="--ex-tab-ref ${prefixes_ref}"
  fi

  # check for computation of reference exon lists
  if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    local with_refs_exons_flag="--with-refs-exons"
  fi

  # check for a defined figure prefix
  local fig_prefix='test'
  if [[ "x${@}" == x*"--fig-prefix"* ]]; then
    local fig_prefix="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local sh_arg=''
  if [[ "x${@}" == x*"--scaleHarm"* ]]; then
    local sh_arg="--scaleHarm"
  fi

  local gene_concat_arg=''
  if [[ "x${@}" == x*"--gene-concat"* ]]; then
    local gene_concat_arg="--gene-concat"
  fi

  local coverage_summary_bank_sel_arg=''
  flag='--cov-bank-sel'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  source ${CONF_FILE}

  local various_length=''
  flag='--various-length'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    vlnb=100
    if [[ -n "$various_length_nbins" ]]; then
      vlnb=${various_length_nbins}
    fi
    local various_length="--various-length $vlnb"

    local side_cut=''
    flag='--side-cut'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local side_cut="--side-cut $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi

    local extensions_arg=''
    flag='--ext-up'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local extensions_arg="${extensions_arg} --ext-up $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi
    flag='--ext-dw'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local extensions_arg="${extensions_arg} --ext-dw $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi

  fi

  if [[ -z "${various_length}" ]]; then
    unset off_set_array; declare -A off_set_array
    off_set_array['3SS']=${cov_out_off}
    off_set_array['5SS']=$(($cov_in_exon-1))
    off_set_array['center']=${cov_center_off}

    unset bed_dir_array; declare -A bed_dir_array
    bed_dir_array['3SS']=${bed_dir}/inEx${cov_in_exon}_outOff${cov_out_off}
    bed_dir_array['5SS']=${bed_dir}/inEx${cov_in_exon}_outOff${cov_out_off}
    bed_dir_array['center']=${bed_dir}/aroundCenter${cov_center_off}

    unset fig_dir_array; declare -A fig_dir_array
    fig_dir_array['3SS']=${fig_dir}/coverage_summary/inEx${cov_in_exon}_outOff${cov_out_off}
    fig_dir_array['5SS']=${fig_dir}/coverage_summary/inEx${cov_in_exon}_outOff${cov_out_off}
    fig_dir_array['center']=${fig_dir}/coverage_summary/aroundCenter${cov_center_off}

    for SS in '3SS' '5SS' 'center'; do
      local bed6_files=$( echo $( echo $( echo ${prefixes_ref} ${input_prefixes} | tr ' ' ',' ) | tr ',' '\n' | awk -v bed_dir=${bed_dir_array[${SS}]}/bed6 -v ss=${SS} '{ print bed_dir"/"$0"_"ss".bed" }' ) | tr ' ' ',' )

      cmd="
      coverage_summary ${bed6_files} ${fig_prefix}_${SS} ${off_set_array[${SS}]} ${fig_dir_array[${SS}]} ${data_dir} ${input_prefixes} ${nbr_clusters} ${min_max_signals} \
        ${prefixes_ref_arg} ${with_refs_exons_flag} ${coverage_summary_bank_sel_arg} ${sh_arg}
      "
      echo "$cmd"
      eval "$cmd"
    done
  else
    local bed6_files=$( echo ${prefixes_ref} ${input_prefixes} | tr ' ' ',' | sed s%'[^,]*'%"${bed_dir}/bed6/&_list.bed"%g )
    local offset=0
    local cov_fig_dir=${fig_dir}/coverage_summary/varLen${vlnb}


    cmd="
    coverage_summary ${bed6_files} ${fig_prefix} ${offset} ${cov_fig_dir} ${data_dir} ${input_prefixes} ${nbr_clusters} ${min_max_signals} \
      ${prefixes_ref_arg} ${with_refs_exons_flag} ${coverage_summary_bank_sel_arg} ${sh_arg} ${various_length} ${side_cut} ${gene_concat_arg} ${extensions_arg}
    "
    echo "$cmd"
    eval "$cmd"
  fi
}
