#!/bin/bash

res_dir=summary_mean_cov_postProc/results
ratio_dir=${res_dir}/ratio
hm_dir=${res_dir}/hm
bxplt_dir=${res_dir}/bxplt
mkdir -vp ${res_dir} ${ratio_dir} ${hm_dir} ${bxplt_dir} >&2

## gather statistics over all the samples
summary_file_all=${ratio_dir}/summary_mean_cov_res.tsv
bash summary_mean_cov_postProc/summary_mean_cov_gatherer.sh ./res > ${summary_file_all}_tmp && \
bash summary_mean_cov_postProc/summary_mean_cov_gatherer.sh ./res_hg38 | tail -n +2 >> ${summary_file_all}_tmp && \
mv ${summary_file_all}_tmp ${summary_file_all}

## select lines for any/byCell results, + remove lines without summary results
summary_file_anyCell=${ratio_dir}/summary_mean_cov_res_anyCell.tsv
head -n 1 ${summary_file_all} > ${summary_file_anyCell}
less ${summary_file_all} | tail -n +2 | grep -vP 'geneVisu_(293T|K562|HeLa|HepG2|MCF7)' | grep -v "NA's" | awk '{ if ( NF == "11" ) print }' >> ${summary_file_anyCell}
summary_file_byCell=${ratio_dir}/summary_mean_cov_res_byCell.tsv
head -n 1 ${summary_file_all} > ${summary_file_byCell}
less ${summary_file_all} | tail -n +2 | grep -P 'geneVisu_(293T|K562|HeLa|HepG2|MCF7)' | grep -v "NA's" | awk '{ if ( NF == "11" ) print }' >> ${summary_file_byCell}


####

ref_order_tab_list_fc() {
  ls report_cov_lists/ChromatinLore_byMark/ChromatinLore_JBC_SF_hgxx.20181210_*.txt
}

marq_str() {
  name="$1"
  echo "$name" | awk -F '[_.]' '{ print $(NF-1) }'
}



#### Exon from any Cell line
signature=anyCell

## select exons only
summary_file=${ratio_dir}/$( basename ${summary_file_anyCell} .tsv )_exon.tsv
head -n 1 ${summary_file_anyCell} > ${summary_file}
less ${summary_file_anyCell} | tail -n +2 | grep -P '^([^\t]*\t){6}[^\t]*_exon\t' >> ${summary_file}

## compute ratios relative to CCE
ratio_tab=${ratio_dir}/$( basename ${summary_file} _exon.tsv )_ratio_exon.tsv
ratio_list='GC_rich_.*exon/AT_rich_.*exon'
ratio_type_list='--ratio-type dnm'
python3 ./summary_mean_cov_postProc/sum_cov_ratioCalc.py ${summary_file} ${ratio_list} ${ratio_type_list} > ${ratio_tab} 2> /dev/null

####
for measure in 'Mean' 'Median'; do
  echo ">>> measure: $measure" >&2
  hm_m_dir=${hm_dir}/${measure}
  bxplt_m_dir=${bxplt_dir}/${measure}
  mkdir -vp ${hm_m_dir} ${bxplt_m_dir} >&2

  ref_order_tab_list="$( ref_order_tab_list_fc )"

  ## build boxplot of ratios
  Rscript ./summary_mean_cov_postProc/ratioTab_2_bxplt.r ${ratio_tab} $( echo $ref_order_tab_list | tr ' ' ',' ) $( echo $( echo "$ref_order_tab_list" | awk -F '[_.]' '{ print $(NF-1) }' | tr '\n' ' ' ) | tr ' ' ',' ) ${bxplt_m_dir}/GC-AT_${signature}_ChromLore --measure ${measure}

  ## build heatmaps of ratios
  for ref_order_tab in $ref_order_tab_list; do
    echo ">>> ref_order_tab: $ref_order_tab" >&2

    marq="$( marq_str ${ref_order_tab} )"
    echo ">>> mark: $marq" >&2

    ratio_selOrder_tab=${ratio_dir}/$( basename $ratio_tab .tsv )_${measure}_selOrder_${marq}.tsv
    python3 ./summary_mean_cov_postProc/selAndOrder.py ${ref_order_tab} ${ratio_tab} > ${ratio_selOrder_tab}
    if [[ "$( wc -l ${ratio_selOrder_tab} )" == 1 ]]; then
      echo '!!! No ratio in '"${ratio_selOrder_tab}" >&2
      continue
    fi

    ratio_hm=${hm_m_dir}/GC-AT_${signature}_ChromLore_${marq}
    echo ">>> heatmap of ratio: $ratio_hm" >&2
    Rscript ./summary_mean_cov_postProc/ratioTab_2_hm.r ${ratio_selOrder_tab} ${ratio_hm} --measure $measure
  done
done


#### Exon by Cell line
signature=byCell

## select exons only
summary_file=${ratio_dir}/$( basename ${summary_file_byCell} .tsv )_exon.tsv
head -n 1 ${summary_file_byCell} > ${summary_file}
less ${summary_file_byCell} | tail -n +2 | grep -P '^([^\t]*\t){6}[^\t]*_exon\t' >> ${summary_file}

## compute ratios relative to CCE
ratio_tab=${ratio_dir}/$( basename ${summary_file} _exon.tsv )_ratio_exon.tsv
ratio_list='GC_rich_.*exon/AT_rich_.*exon'
ratio_type_list='--ratio-type dnm'
python3 ./summary_mean_cov_postProc/sum_cov_ratioCalc.py ${summary_file} ${ratio_list} ${ratio_type_list} > ${ratio_tab}

####
for measure in 'Mean' 'Median'; do
  echo ">>> measure: $measure" >&2
  hm_m_dir=${hm_dir}/${measure}
  bxplt_m_dir=${bxplt_dir}/${measure}
  mkdir -vp ${hm_m_dir} ${bxplt_m_dir} >&2

  ref_order_tab_list="$( ref_order_tab_list_fc )"

  ## build boxplot of ratios
  Rscript ./summary_mean_cov_postProc/ratioTab_2_bxplt.r ${ratio_tab} $( echo $ref_order_tab_list | tr ' ' ',' ) $( echo $( echo "$ref_order_tab_list" | awk -F '[_.]' '{ print $(NF-1) }' | tr '\n' ' ' ) | tr ' ' ',' ) ${bxplt_m_dir}/GC-AT_${signature}_ChromLore --measure ${measure}

  ## build heatmaps of ratios
  for ref_order_tab in $ref_order_tab_list; do
    echo ">>> ref_order_tab: $ref_order_tab" >&2

    marq="$( marq_str ${ref_order_tab} )"
    echo ">>> mark: $marq" >&2

    ratio_selOrder_tab=${ratio_dir}/$( basename $ratio_tab .tsv )_${measure}_selOrder_${marq}.tsv
    python3 ./summary_mean_cov_postProc/selAndOrder.py ${ref_order_tab} ${ratio_tab} > ${ratio_selOrder_tab}
    if [[ "$( wc -l ${ratio_selOrder_tab} )" == 1 ]]; then
      echo '!!! No ratio in '"${ratio_selOrder_tab}" >&2
      continue
    fi

    ratio_hm=${hm_m_dir}/GC-AT_${signature}_ChromLore_${marq}
    echo ">>> heatmap of ratio: $ratio_hm" >&2
    Rscript ./summary_mean_cov_postProc/ratioTab_2_hm.r ${ratio_selOrder_tab} ${ratio_hm} --measure $measure
  done
done
