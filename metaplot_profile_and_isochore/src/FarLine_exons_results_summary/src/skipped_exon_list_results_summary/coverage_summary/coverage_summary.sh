#!/bin/bash

coverage_summary_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${coverage_summary_dir}/metagene_coverage.sh
source ${coverage_summary_dir}/metagene_diffCoverage.sh
source ${coverage_summary_dir}/kmeans_clustering.sh
source ${coverage_summary_dir}/coverage_summary_stats.sh

coverage_summary () {
  ## process arguments
  local bed6_files=$1 #~ './results/splicing_site_regions/inEx200_outOff500/bed6/exon_up_siDNMT3b-siGL2_3SS.bed'
  local sign=$2 #~ '3SS'
  local off_set=$3 #~ 500
  local out_dir=$4 #~ ./results/figures/coverage_summary
  local data_dir=$5
  local prefixes_query=$6
  local nbr_clusters=$7
  local min_max_signals=$8
  # local prefixes_ref=$9

  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local prefixes_ref=''
  flag='--ex-tab-ref'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local prefixes_ref="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  local sh_arg=''
  if [[ "x${@}" == x*"--scaleHarm"* ]]; then
    local sh_arg="--scaleHarm"
  fi

  local various_length=''
  flag='--various-length'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local various_length="--various-length $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"

    local side_cut=''
    flag='--side-cut'
    if [[ "x${@}" == x*"${flag}"* ]]; then
      local side_cut="--side-cut $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
    fi

    local extensions_arg=''
    flag='--ext-up'
    if [[ -n "${ext_up}" ]]; then
      local extensions_arg="${extensions_arg} --ext-up ${ext_up}"
    fi
    flag='--ext-dw'
    if [[ -n "${ext_dw}" ]]; then
      local extensions_arg="${extensions_arg} --ext-dw ${ext_dw}"
    fi
  fi

  local gene_concat_arg=''
  if [[ "x${@}" == x*"--gene-concat"* ]]; then
    local gene_concat_arg="--gene-concat"
  fi

  # check for computation of ASE and CE lists
  bed6_file_array=( $( echo ${bed6_files} | tr ',' ' ' ) )
  if [[ "x${@}" != x*"--with-refs-exons"* ]]; then
    bed6_file_text="$( echo ${bed6_file_array[*]} | tr ' ' '\n' )"
    for prefix_ref in $( echo ${prefixes_ref} | tr ',' ' ' ); do
      bed6_file_text="$( echo "${bed6_file_text}" | grep -v ${prefix_ref} )"
    done
    bed6_file_array=( $( echo ${bed6_file_text} | tr '\n' ' ' ) )
  fi


  # gather the prefixes
  local prefixes=$( echo ${prefixes_ref} ${prefixes_query} | tr ' ' ',' )

  ## get list of selected coverage tracks if some specified
  # coverage_summary_bank_sel=/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/cov_sel.txt
  coverage_summary_bank_sel=''
  flag='--cov-bank-sel'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    local coverage_summary_bank_sel="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  coverage_summary_bank_dir=${coverage_summary_dir}/coverage_bank
  for expSeq_source in $( find ${coverage_summary_bank_dir}/ -name "*_script.sh" | sort ); do
    # echo -e "\t!!! $( basename ${expSeq_source} )"
    ## check if coverage track is among selected tracks
    if [[ -n "$coverage_summary_bank_sel" ]]; then
      # echo "--- $coverage_summary_bank_sel" >&2
      check_sel_script="$( cat $( echo ${coverage_summary_bank_sel} | tr ',' ' ' ) | grep -vE '^#' | grep -P "^$( basename $expSeq_source )" )"
      # echo -e "\t>>> $( basename $expSeq_source )" >&2
      if [[ -z "$check_sel_script" ]]; then
        continue
      fi
    fi

    echo -e "\t>>> $( basename ${expSeq_source} )"
    mean_ymin='-' mean_ymax='-' median_ymin='-' median_ymax='-'
    source ${expSeq_source}

    ## get y-limits if specified in file of selected coverage tracks
    col_lim_mean="$( echo "$check_sel_script" | awk '{ print $2 }' )"
    if [[ "x$col_lim_mean" == x*","* ]]; then
      mean_ymin="$( echo "$col_lim_mean" | cut -d ',' -f 1 )"
      mean_ymax="$( echo "$col_lim_mean" | cut -d ',' -f 2 )"
    fi
    col_lim_median="$( echo "$check_sel_script" | awk '{ print $3 }' )"
    if [[ "x$col_lim_median" == x*","* ]]; then
      median_ymin="$( echo "$col_lim_median" | cut -d ',' -f 1 )"
      median_ymax="$( echo "$col_lim_median" | cut -d ',' -f 2 )"
    fi

    # local ref_chr_arg=''
    # local sum_cov_chr_arg=''
    # if [[ "x${check_sel_script}" == x*'--ref-add-chr'* ]]; then
    #   local ref_chr_arg='--ref-add-chr'
    #   local sum_cov_chr_arg='--add-chr 1'
    # fi


    # check for two condition are present for a differential coverage computing
    local nbr_comp_pair=$( echo ${comp_pair} | awk -F ',' '{ print NF }' )
    local kmean_sig_type=$( echo $comp_pair | cut -d ',' -f 1 )
    local diffCoverage_bool=0
    if [[ "x${nbr_comp_pair}" == x"2" ]]; then
      diffCoverage_bool=1
      kmean_sig_type="${kmean_sig_type} ${comp_pair}"
    fi

    local expSeq_sign=$( basename ${expSeq_source} _script.sh )
    # echo "!!! $expSeq_source"
    # echo "!!! $expSeq_sign"

    ## test if human refs with 'chr'
    bw_info="$( Rscript ${coverage_summary_dir}/checkBigWig.r $( echo ${expSeq_bw_files} | cut -d ',' -f 1 ) )"
    bw_info="$( echo "$bw_info" | grep -e chrom -A 1 | grep -e chr[0-9XYM] )"
    echo "--- ${expSeq_bw_files}" >&2
    echo "--- ${bw_info}" >&2
    if [ -z "$bw_info" ] ; then
    echo "--- no chr" >&2
      echo "ref as X"
      ref_chr_arg=""
    else
    echo "--- with chr" >&2
      echo "ref as chrX"
      ref_chr_arg="--ref-add-chr"
    fi

    ## check enough slots for submitting a qsub job
    nbr_jobs="$( qstat | wc -l )"
    while [[ "$nbr_jobs" -gt 14500 ]]; do
      sleep 60
      nbr_jobs="$( qstat | wc -l )"
    done


    echo "metagene profile" >&2
    echo "coverage" >&2
    local expSeq_out_dir=${out_dir}/${expSeq_sign}/cov
    mkdir -vp ${expSeq_out_dir}
    cmd="metagene_coverage ${expSeq_bw_files} ${expSeq_cond} ${expSeq_rep} ${bed6_files} ${sign} ${expSeq_out_dir} ${off_set} $( echo $comp_pair | cut -d ',' -f 1 ) ${prefixes} \
      ${color_pallette} --hm-ymax ${hm_ymax} --mean-ymax ${mean_ymin},${mean_ymax} --median-ymax ${median_ymin},${median_ymax} ${sh_arg} ${ref_chr_arg} ${various_length} ${side_cut} ${gene_concat_arg} ${extensions_arg}"
    echo "$cmd" | tr ' ' '\n'
    eval "$cmd"

    # compute mean signal per annotation
    cmd="source ${coverage_summary_dir}/coverage_means.sh; coverage_means ${expSeq_bw_files} ${bed6_files} ${prefixes} ${expSeq_out_dir} ${sign} ${color_pallette} ${ref_chr_arg} ${various_length} ${side_cut} ${gene_concat_arg} ${extensions_arg}"
    qsub_script=${expSeq_out_dir}/qsub_script_${sign}_means.sh
    less ${coverage_summary_dir}/run_cov_means.sh | sed s%'ppp_cmd_ppp'%"$cmd"% | sed s%'ppp_logs_[erout]*_dir_ppp'%"${expSeq_out_dir}"% > $qsub_script
    cmd="qsub $qsub_script"
    echo "$cmd" | tr ' ' '\n'
    eval "$cmd"

    if [[ -n "$various_length" ]]; then
      out_summary_cov=${expSeq_out_dir}/summary_mean_cov.tsv
      cmd="source ${coverage_summary_dir}/coverage_summary_stats.sh; coverage_summary_stats ${expSeq_bw_files} ${expSeq_cond} ${expSeq_rep} ${bed6_files} ${prefixes} ${out_summary_cov} ${ref_chr_arg}"

      qsub_script=${expSeq_out_dir}/qsub_script_stat.sh
      less ${coverage_summary_dir}/run_cov_sum_stats.sh | sed s%'ppp_cmd_ppp'%"$cmd"% | sed s%'ppp_logs_[erout]*_dir_ppp'%"${expSeq_out_dir}"% > $qsub_script
      cmd="qsub $qsub_script"
      echo "$cmd"
      eval "$cmd"
    fi


    if [[ "x${diffCoverage_bool}" == x"1" ]]; then
      echo "diffCoverage" >&2
      local expSeq_out_dir=${out_dir}/${expSeq_sign}/diffCov
      mkdir -vp ${expSeq_out_dir}
      # metagene_diffCoverage ${expSeq_bw_files} ${expSeq_cond} ${expSeq_rep} ${bed6_files} ${sign} ${expSeq_out_dir} ${off_set} ${comp_pair} ${prefixes} \
        # ${color_pallette} ${sh_arg}
    fi

    # echo "K-means clustering" >&2
    # for nbr_cluster in $( echo ${nbr_clusters} | tr ',' ' ' ); do
    #   for min_max_signal in $( echo ${min_max_signals} | tr ',' ' ' ); do
    #     unset expSeq_out_dir_array; declare -A expSeq_out_dir_array
    #     expSeq_out_dir_array[$( echo $comp_pair | cut -d ',' -f 1 )]=${out_dir}/${expSeq_sign}/cov/${nbr_cluster}clusters/${min_max_signal}mms
    #     if [[ "x${nbr_comp_pair}" == x"2" ]]; then
    #       expSeq_out_dir_array[$comp_pair]=${out_dir}/${expSeq_sign}/diffCov/${nbr_cluster}clusters/${min_max_signal}mms
    #     fi
    #
    #     for bed6_file in ${bed6_file_array[*]}; do
    #       # echo "--- bed6_file: $bed6_file" >&2
    #       for sig_type in ${kmean_sig_type}; do
    #         # echo "--- working signal: ${sig_type}" >&2
    #         mkdir -vp ${expSeq_out_dir_array[$sig_type]}
    #         kmeans_clustering ${expSeq_bw_files} ${expSeq_cond} ${expSeq_rep} ${bed6_file} $( basename ${bed6_file} .bed ) ${expSeq_out_dir_array[$sig_type]} ${off_set} ${sig_type} ${nbr_cluster} ${min_max_signal} --hm-ymax ${ymax}
    #       done
    #     done
    #   done
    # done

  done

}
