#!/bin/bash



res_dir=summary_mean_cov_postProc/results
ratio_dir=${res_dir}/ratio
hm_dir=${res_dir}/hm
mkdir -vp ${res_dir} ${ratio_dir} ${hm_dir} >&2

## gather statistics over all the samples
summary_file_all=${ratio_dir}/summary_mean_cov_res.tsv
bash summary_mean_cov_postProc/summary_mean_cov_gatherer.sh ./res > ${summary_file_all}
bash summary_mean_cov_postProc/summary_mean_cov_gatherer.sh ./res_hg38 | tail -n +2 >> ${summary_file_all}

## select exons only
summary_file=${ratio_dir}/summary_mean_cov_res_exon.tsv
less ${summary_file_all} | grep -P '^(([^\t]*\t){6}[^\t]*_exon\t|Min)' > ${summary_file}


####
ref_order_tab_list_fc() {
  ls ./report_cov_lists/{H,cistrome_H,histone,cistrome_histone}* #| grep -P '^\./report_cov_lists/(cistrome_histone|histone)'
}

marq_str() {
  local ref_order_tab=$1
  if [[ "$( basename ${ref_order_tab} )" == "H"?_* || "$( basename ${ref_order_tab} )" == "cistrome_H"* || "$( basename ${ref_order_tab} )" == "cistrome_histone"* || "$( basename ${ref_order_tab} )" == "histone"* ]]; then
    marq=$( basename ${ref_order_tab} )
  else
    marq=$( basename ${ref_order_tab} | cut -d '_' -f 1 )
  fi

  echo "$marq"
}

for measure in 'Mean' 'Median'; do
  hm_m_dir=${hm_dir}/${measure}
  mkdir -vp ${hm_m_dir} >&2

  ####
  ## compute ratios relative to CCE
  ratio_tab=${ratio_dir}/summary_mean_cov_ratioToCCE_exon.tsv
  ratio_list='GC_rich_exon/CCE_exon,AT_rich_exon/CCE_exon'
  ratio_type_list='--ratio-type dn,dn'
  python3 ./summary_mean_cov_postProc/sum_cov_ratioCalc.py ${summary_file} ${ratio_list} ${ratio_type_list} > ${ratio_tab}

  ## build heatmaps of ratios
  for ref_order_tab in $( ref_order_tab_list_fc ); do
    echo "${ref_order_tab}" >&2
    # ref_order_tab=./report_cov_lists/histoneMeth_complete_refined
    marq="$( marq_str ${ref_order_tab} )"
    echo "$marq"

    ratio_selOrder_tab=${ratio_dir}/$( basename $ratio_tab .tsv )_selOrder_${marq}.tsv
    python3 ./summary_mean_cov_postProc/selAndOrder.py ${ref_order_tab} ${ratio_tab} > ${ratio_selOrder_tab}

    ratio_hm=${hm_m_dir}/ratiosToCCE_exon_${marq}.png
    echo $ratio_hm
    Rscript ./summary_mean_cov_postProc/ratioTab_2_hm.r ${ratio_selOrder_tab} ${ratio_hm} --measure $measure
  done

  ####
  ## compute ratios GC-rich exons over AT-rich exons
  ratio_tab=${ratio_dir}/summary_mean_cov_ratioGCvsAT_exon.tsv
  ratio_list='GC_rich_exon/AT_rich_exon'
  ratio_type_list='--ratio-type dnm'
  python3 ./summary_mean_cov_postProc/sum_cov_ratioCalc.py ${summary_file} ${ratio_list} ${ratio_type_list} > ${ratio_tab}

  ## build heatmaps of ratios
  for ref_order_tab in $( ref_order_tab_list_fc ); do
    echo "${ref_order_tab}" >&2
    # ref_order_tab=./report_cov_lists/histoneMeth_complete_refined
    marq="$( marq_str ${ref_order_tab} )"
    echo "$marq"

    ratio_selOrder_tab=${ratio_dir}/$( basename $ratio_tab .tsv )_selOrder_${marq}.tsv
    python3 ./summary_mean_cov_postProc/selAndOrder.py ${ref_order_tab} ${ratio_tab} > ${ratio_selOrder_tab}

    ratio_hm=${hm_m_dir}/ratiosGCvsAT_exon_${marq}.png
    echo $ratio_hm
    Rscript ./summary_mean_cov_postProc/ratioTab_2_hm.r ${ratio_selOrder_tab} ${ratio_hm} --measure $measure
  done
done
