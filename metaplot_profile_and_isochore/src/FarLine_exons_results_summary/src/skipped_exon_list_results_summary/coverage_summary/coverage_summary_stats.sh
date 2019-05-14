#!/bin/bash

coverage_summary_stats_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

coverage_summary_stats () {
  ## process arguments
  local expSeq_bw_files=$1
  local expSeq_cond=$2 
  local expSeq_rep=$3 
  local bed6_files=$4 
  local prefixes=$5
  local out_summary_cov=$6

  local sum_cov_chr_arg=''
  flag='--ref-add-chr'
  if [[ "x${@}" == x*"${flag}"* ]]; then
      local sum_cov_chr_arg='--add-chr 1'
  fi



  ebf_array=( $( echo ${expSeq_bw_files} | tr ',' ' ' ) )
  econd_array=( $( echo ${expSeq_cond} | tr ',' ' ' ) )
  erep_array=( $( echo ${expSeq_rep} | tr ',' ' ' ) )

  bed6_array=( $( echo ${bed6_files} | tr ',' ' ' ) )
  prefixes_array=( $( echo ${prefixes} | tr ',' ' ' ) )

  echo -e "Min\t1st_Qu\tMedian\tMean\t3rd_Qu\tMax\tannot_list\tcond\trep" > ${out_summary_cov}

  for bed_idx in $( seq 1 ${#bed6_array[*]} ); do
    real_bed_idx=$((bed_idx-1))
    for cov_idx in $( seq 1 ${#ebf_array[*]} ); do
      real_cov_idx=$((cov_idx-1))
      cmd="python3 ${coverage_summary_stats_dir}/metagene_coverage_functions/retrieve_bw_values_mean.py ${ebf_array[$real_cov_idx]} ${bed6_array[$real_bed_idx]} - --stdout ${sum_cov_chr_arg} | Rscript ${coverage_summary_stats_dir}/metagene_coverage_functions/summary_stats.r"
      echo "$cmd"
      sum_res="$( eval "$cmd" | tail -n +2 | sed s%'\(^\ *\|\ *$\)'%''%g | sed s%'\ \ *'%'\t'%g )"
      echo "${sum_res}" ${prefixes_array[$real_bed_idx]} ${econd_array[$real_cov_idx]} ${erep_array[$real_cov_idx]} | tr ' ' '\t' >> ${out_summary_cov}
    done
  done
}
