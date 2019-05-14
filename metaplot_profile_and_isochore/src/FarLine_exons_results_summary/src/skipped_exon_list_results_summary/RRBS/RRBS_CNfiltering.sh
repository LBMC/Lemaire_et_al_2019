#!/bin/bash

RRBS_CNfiltering_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cn2distFile () {
  RRBS_cn_dir=$2
  cn_distToSS_dir=$3
  SS_bed_base=$4
  echo ${1} | sed s,${RRBS_cn_dir},${cn_distToSS_dir},g | sed s/.cn/_${SS_bed_base}_distToSS.cn/g
}

RRBS_CNfiltering () {
  ## process args
  local prefixes=$1 #~ 'exon_up_siDNMT3b-siGL2'
  local SS=$2 #~ '3SS'
  local exons_bed_dir=$3 #~ ${out_dir}/exons_lists/bed6
  local ss_bed_dir=$4
  local rate_limits=$5 #~ 13,87
  local RRBS_cn_files=$6
  local cn_distToSS_dir=$7
  local fract_out_dir=$8
  local start=$9
  local end=${10}
  local window_size=${11}
  local RRBS_cond=${12}
  local RRBS_rep=${13}

  local NBRCORES=1
  if [[ "x$@" == x*"--nbr-cores"* ]]; then
    NBRCORES=$( echo "$@" | awk -F '--nbr-cores' '{ print $NF }' | cut -d ' ' -f 2 )
  fi


  ## filtered the CN files for the regions of interest
  for prefix in $( echo ${prefixes} | tr ',' ' ' ); do
    local exons_bed="${exons_bed_dir}/${prefix}_list.bed"
    local SS_bed_base=${prefix}_${SS}
    local SS_bed="${ss_bed_dir}/${SS_bed_base}_ext.bed"

    local RRBS_distToSS_files=''
    local RRBS_cn_files_array=( $( echo ${RRBS_cn_files} | tr ',' ' ' ) )

    local nbr_proc=0
    local nbr_task=0
    for cn_input in ${RRBS_cn_files_array[*]}; do
      # filter for the cytosine in the splicing regions, and compute distance from the cytosine to the splicing site
      local RRBS_cn_dir=$( dirname ${cn_input} )
      local cn_distToSS_file=$( cn2distFile ${cn_input} ${RRBS_cn_dir} ${cn_distToSS_dir} ${SS_bed_base} ).gz #~ ${cn_distToSS_dir}/$( basename ${cn_input} .cn )_${SS_bed_base}_distToSS.cn
      RRBS_distToSS_files=${RRBS_distToSS_files},${cn_distToSS_file}

      echo ">>> filtering and distToSS" >&2
      python3 ${RRBS_CNfiltering_dir}/CNIntersectBed.py ${cn_input} ${SS_bed} --annot-colname 'exon_id' --annot-id 4 2> /dev/null | \
        awk '{OFS="\t"; if ($NF!="NA") print}' | \
        python3 ${RRBS_CNfiltering_dir}/groupToLines.py | \
      Rscript ${RRBS_CNfiltering_dir}/mCPosTodistToSS.r - ${exons_bed} ${SS} | \
      gzip -c > ${cn_distToSS_file} &

      nbr_proc=$(($nbr_proc+1))
      nbr_task=$(($nbr_task+1))
      if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#RRBS_cn_files_array[@]}" ]]; then
        nbr_proc=0
        wait
      fi


    done
    RRBS_distToSS_files=$( echo ${RRBS_distToSS_files} | sed s/^','// )

    # compute the fractions per 5mC rate categories
    echo ">>> fraction" >&2
    local fract_out_file=${fract_out_dir}/${SS_bed_base}_fract.RDS
    Rscript ${RRBS_CNfiltering_dir}/cnDist2fraction.r ${RRBS_distToSS_files} ${rate_limits} -${start},${end} ${window_size} ${fract_out_file} ${RRBS_cond} ${RRBS_rep}
    # exit
  done

}
