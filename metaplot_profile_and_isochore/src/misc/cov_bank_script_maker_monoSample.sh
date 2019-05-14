

data_dir=/home/slemaire/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data
bw_dir=$1 #/home/slemaire/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data/bw_files/BIGWIG_hg38
sub_bw_dir="$( echo $bw_dir | sed s%"$data_dir/"%""% )"
cb_dir=/home/slemaire/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/coverage_bank


script_maker () {

  echo '#!/bin/bash'

  sign=$1
  sub_bw_dir=$2

  echo 'public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  expSeq_bw_dir=${data_dir}/'${sub_bw_dir}'
  encode_sign='${sign}'
  source ${public_script_dir}/public_common_encode.sh
  '
}


for xxx in $( ls $bw_dir/*.bw ); do
  bw_name=$( basename $xxx .bw )
  # echo $bw_name

  script=${cb_dir}/${bw_name}_script.sh
  basename $script
  script_maker $bw_name $sub_bw_dir > $script

  # break
done


####
