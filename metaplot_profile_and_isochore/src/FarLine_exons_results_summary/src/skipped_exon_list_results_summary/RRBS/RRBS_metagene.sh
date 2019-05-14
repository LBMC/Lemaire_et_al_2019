#!/bin/bash

RRBS_metagene_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cn2distFile () {
  RRBS_cn_dir=$2
  cn_distToSS_dir=$3
  SS_bed_base=$4
  echo ${1} | sed s,${RRBS_cn_dir},${cn_distToSS_dir},g | sed s/.cn/_${SS_bed_base}_distToSS.cn/g
}

RRBS_metagene () {
  ## process args
  local prefixes=$1 #~ 'exon_up_siDNMT3b-siGL2'
  local SS=$2 #~ '3SS'
  local out_fig_prefix=$3 #~ ./results/figures/inEx200_outOff500/RRBS/methCyt/win10b/exon_up_siDNMT3b-siGL2_3SS
  local fract_out_dir=$4
  local comp_pair=$5
  local mode=$6

  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  ## build the metagenes of RRBS (median on the replicates)
  fract_rds_list=$( echo $( echo ${prefixes} | tr ',' '\n' | sed s,.*,"${fract_out_dir}/&_${SS}_fract.RDS", ) | tr ' ' ',' )
  Rscript ${RRBS_metagene_dir}/RRBS_metagene.r ${fract_rds_list} ${out_fig_prefix} --prefixes ${prefixes} --comp-pair ${comp_pair} --mode ${mode} ${color_pallette}

}
