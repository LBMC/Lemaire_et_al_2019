#!/bin/bash

nb_ref_lists_fc () {
  local exon_ref_list=$1
  local exon_ref_nb_list=$2
  local nb_type=$3

  exon_ref_array=( $( echo "$exon_ref_list" | tr ',' ' ' ) )
  exon_ref_nb_array=( $( echo "$exon_ref_nb_list" | tr ',' ' ' ) )

  if [[ "x${nb_type}" == x"nbup" ]]; then
    trans_val="-1"
  elif [[ "x${nb_type}" == x"nbdown" ]]; then
    trans_val="1"
  fi

  for (( iii=0; iii<${#exon_ref_array[@]}; iii++ )); do
    exon_ref=${exon_ref_array[$iii]}
    exon_ref_nb=${exon_ref_nb_array[$iii]}

    head -1 ${exon_ref} > ${exon_ref_nb}
    tail -n +2 ${exon_ref} | awk -F '_' -v trans_val=${trans_val} '{ print $1"_"trans_val+$2 }' >> ${exon_ref_nb}
  done
}
