#!/bin/bash

rna_folding_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rna_folding () {
  local fasta_file_list=$1
  local names=$2
  local out_rna_folding_dir=$3
  local out_fig_dir=$4

  mkdir -vp "${out_rna_folding_dir}" "${out_fig_dir}" >&2

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

  ## compute the foldings and the associated Minimum Free Energies (MFEs).
  unset res_file_array; declare res_file_array
  for fasta_file in $( echo "${fasta_file_list}" | tr ',' ' ' ); do
    local res_file=${out_rna_folding_dir}/$( basename ${fasta_file} .fa ).txt
    res_file_array=( ${res_file_array[*]} ${res_file} )

    cmd="
    ${rna_folding_dir}/ViennaRNA-2.4.4/bin/RNAfold --noPS ${fasta_file} > ${res_file}
    "
    echo "$cmd"
    eval "$cmd" &

    nbr_proc=$(($nbr_proc+1))
    nbr_task=$(($nbr_task+1))
    if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
      nbr_proc=0
      wait
    fi

  done
  wait

  ## make repartition plots of MFE
  cmd="
  Rscript ${rna_folding_dir}/MFE_repart.r \
    $( echo "${res_file_array[*]}" | tr ' ' ',' ) \
    ${names} \
    ${out_fig_dir} \
    ${fig_prefix} \
    ${color_pallette}
  "
  echo "$cmd"
  eval "$cmd"
}
