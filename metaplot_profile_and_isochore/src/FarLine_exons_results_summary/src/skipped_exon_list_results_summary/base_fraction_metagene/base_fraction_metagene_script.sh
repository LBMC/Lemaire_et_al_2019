#!/bin/bash

base_fraction_metagene_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${base_fraction_metagene_script_dir}/base_fraction_metagene.sh

base_fraction_metagene_script () {
  local exon_table_prefixes=$1
  local FASTA_DIR=$2
  local res_dir=$3
  local CONF_FILE=$4
  # local exon_table_prefixes_ref=$5

  local exon_table_prefixes_ref=''
  flag='--ex-tab-ref'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local exon_table_prefixes_ref="--ex-tab-ref $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # get subparts of base fraction if specified
  local bf_arg='--bf-parts --bf-all'
  flag='--bf-parts'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local bf_arg="--bf-parts $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  source ${CONF_FILE}

  # check for computation of ref exon lists
  if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    local with_refs_exons_flag="--with-refs-exons"
  fi

  # check for a defined figure prefix
  local fig_prefix='test'
  if [[ "x${@}" == x*"--fig-prefix"* ]]; then
    local fig_prefix="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # check for specific colors in metagene plots
  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi


  # compute base fraction for all the region relative definition, splicing sites, and window size
  in_exon_array=( $( echo ${in_exon_list} | tr ',' ' ' ) )
  out_off_array=( $( echo ${out_off_list} | tr ',' ' ' ) )
  for ((iii=0; iii<${#in_exon_array[@]}; iii++)); do
    in_exon=${in_exon_array[$iii]}
    out_off=${out_off_array[$iii]}
    for SS in '3SS' '5SS'; do
      for window_size in $( echo ${window_sizes} | tr ',' ' ' ); do
        # echo "!!! window_size: ${window_size}b" >&2
        local metagene_param="$window_size,$in_exon,$out_off"
        local bfm_out_dir=${res_dir}/splicing_site_regions/inEx${in_exon}_outOff${out_off}/win${window_size}b

        cmd="
        base_fraction_metagene ${exon_table_prefixes} ${FASTA_DIR} ${SS} ${metagene_param} ${bfm_out_dir} ${min_reads}\
          ${exon_table_prefixes_ref} --nbr-cores ${NBRCORES} ${with_refs_exons_flag} --fig-prefix ${fig_prefix}_${SS} ${color_pallette} ${bf_arg}
        "
        echo "$cmd"
        eval "$cmd"
      done
    done
  done
}
