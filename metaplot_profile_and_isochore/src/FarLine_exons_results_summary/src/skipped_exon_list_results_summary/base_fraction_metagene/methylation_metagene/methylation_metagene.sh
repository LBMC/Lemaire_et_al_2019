#/bin/bash

methylation_metagene_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

methylation_metagene () {
  # local CONF_FILE=$6
  local window_size=$6
  local in_exon=$7
  local out_off=$8
  local min_reads=${12}

  local fig_dir=$5

  local exon_seq_files=$1
  local meth_files=$2
  local exon_seq_names=$3
  local SS_type=$4 #~ 3SS or 5SS

  local expSeq_cond=$9
  local expSeq_rep=${10}
  local comp_pair=${11}


  local NBRCORES=1
  local nbr_proc=0
  local nbr_task=0
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


  local pos_index=1

  unset start_rel; declare -A start_rel
  local start_rel_script=${methylation_metagene_dir}/start_rel.r
  local start_rel['3SS']=$( Rscript ${start_rel_script} "$(($pos_index-$out_off))" $window_size 1 )
  local start_rel['5SS']=$( Rscript ${start_rel_script} "-$(($pos_index+$in_exon-1))" $window_size 1 )

  local out_prefix=${fig_dir}/${fig_prefix}

  ## determine the running mode to consider
  mode_array=( '--raw-mode' )
  comp_pair_array=( $( echo ${comp_pair} | tr ',' ' ' ) )
  if [ "${#comp_pair_array[@]}" -gt "1" ]; then
    mode_array=( ${mode_array[@]} '--diff-mode' )
  elif [ "${#comp_pair_array[@]}" -gt "1" ]; then
    mode_array=()
    echo "!!! 'comp_pair' is empty." >&2
  fi

  for sig_mode in ${mode_array[*]}; do
    Rscript ${methylation_metagene_dir}/meth_meta_plot.r \
      ${meth_files} \
      ${exon_seq_files} \
      ${exon_seq_names} \
      ${expSeq_cond} \
      ${expSeq_rep} \
      ${comp_pair} \
      ${start_rel[$SS_type]} \
      ${window_size} \
      ${min_reads} \
      ${out_prefix} \
      ${color_pallette} \
      ${sig_mode}
  done
}
