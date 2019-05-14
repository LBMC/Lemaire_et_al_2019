#!/bin/bash

plot2D_deltaPsi_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

plot2D_deltaPsi() {
  all_data_query=$1 #~ /home/sebastien/work_temp/id_card_temp/results/all_data_tab/all_data_siPP_siGL2.tab
  all_data_ref=$2 #~ /home/sebastien/work_temp/id_card_temp/results/all_data_tab/all_data_siDNMT3b-siGL2.tab
  query_name=$3
  ref_name=$4
  out_dir=$5 #~ ./bouh_fig

  mkdir -vp ${out_dir} >&2

  tab_names="${ref_name},${query_name}"


  unset corr_array; declare -A corr_array
  corr_array['all']=''
  corr_array['corr']='--corr'
  corr_array['antiC']='--antiC'

  fig_name_list=''
  for corr_val in 'all' 'corr' 'antiC'; do
    suffix="${query_name}_vs_${ref_name}_${corr_val}"
    local fig_name="$(
      Rscript ${plot2D_deltaPsi_dir}/comp_deltaPSI_siRef_other.r \
      ${all_data_ref} \
      ${all_data_query} \
      ${tab_names} \
      --psi 0.1 \
      --pval 0.05 \
      --suf ${suffix} \
      --out_dir ${out_dir} \
      --no-filt-ref \
      --common-only \
      --pval-highlight \
      --no-corr \
      --no-enrich \
      ${corr_array[$corr_val]} \
      ;
    )"
    fig_name_list=( "${fig_name_list[*]}" "${fig_name}" )
  done
  echo ${fig_name_list[*]} | tr ' ' ','
}
