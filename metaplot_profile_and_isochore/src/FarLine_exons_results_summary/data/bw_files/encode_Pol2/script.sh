

script_maker () {
  sub_bw_dir=$1
  encode_sign=$2

  txt='#!/bin/bash

public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

expSeq_bw_dir=${data_dir}/'${sub_bw_dir}'
encode_sign='${encode_sign}'
source ${public_script_dir}/public_common_encode.sh
'
  echo "$txt"
}



data_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data'
bw_dir=bw_files/encode_Pol2
cov_scr_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/coverage_bank'
csvs_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/cov_sum_bank_sel'
target_list="eGFP-POLR2H
POLR2A
POLR2AphosphoS2
POLR2AphosphoS5
POLR2B
POLR2G
"

for target in $target_list; do
  link_dir=${bw_dir}/${target}
  mkdir -vp $link_dir >&2

  bw_list=$( ls ~/analyses_SEBASTIEN/data/bigwig_files/encode_Pol2/*_${target}_*.bw )
  echo ">>> $target"
  # echo "$bw_list"

  for bw in $bw_list; do
    cmd="ln -s $bw ${link_dir}/$( basename ${bw} )"
    # echo "$cmd"
    # eval "$cmd" 2> /dev/null

    encode_sign="$( basename $bw .bw )"

    out_script=${cov_scr_dir}/${encode_sign}_script.sh
    basename $out_script
    # script_maker ${link_dir} ${encode_sign} > ${out_script}
  done
  # break
done | \
while read line; do
  # if [[ "$line" == ">>> "* ]]; then
  #   out_csvs=${csvs_dir}/$( echo "$line" | awk '{ print $NF }' )_tracks.txt
  #   echo ${out_csvs}
  #   echo '' > ${out_csvs}
  # else
  #   echo "$line" >> ${out_csvs}
  # fi
done | less
