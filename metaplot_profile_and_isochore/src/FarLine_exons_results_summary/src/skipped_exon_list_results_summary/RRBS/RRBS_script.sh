#!/bin/bash

RRBS_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${RRBS_script_dir}/RRBS_CNfiltering.sh
source ${RRBS_script_dir}/RRBS_metagene.sh

RRBS_script () {
  local exon_table_prefixes=$1
  local res_dir=$2
  local data_dir=$3
  local exons_bed_dir=$4 #~ ${exons_list_dir}/bed6
  local ss_dir=$5 #~ ${ss_list_dir}/bed6
  local CONF_FILE=$6
  local fig_dir=$7 #~ ./results/figures
  local exon_table_prefixes_ref=$8

  # check for a defined figure prefix
  local fig_prefix='test'
  if [[ "x${@}" == x*"--fig-prefix"* ]]; then
    local fig_prefix="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # check for defined color for the figure
  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  source ${CONF_FILE}

  ## define some variables
  # directory of input cn files and other experimental details
  local cn_file_prefix=CG_RRBS_MCF7
  local RRBS_cn_dir="${data_dir}/CN_files/RRBS/CG_pos"
  local RRBS_cn_files="$( echo $( find ${RRBS_cn_dir}/ -maxdepth 1 -name "${cn_file_prefix}siGL2n?.cn" | sort ) $( find ${RRBS_cn_dir}/ -maxdepth 1 -name "${cn_file_prefix}siPPn?.cn" | sort ) | tr ' ' ',' )"
  local RRBS_cond=$( echo $( echo ${RRBS_cn_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | sed s/${cn_file_prefix}// | cut -d 'n' -f 1 ) | tr ' ' ',' )
  local RRBS_rep=$( echo $( echo ${RRBS_cn_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | sed s/${cn_file_prefix}// | cut -d 'n' -f 2 | cut -d '.' -f 1 ) | tr ' ' ',' )
  local comp_pair=$( echo ${RRBS_cond} | tr ',' '\n' | uniq | tr '\n' ',' )

  unset start_array; declare -A start_array
  start_array['3SS']=${cov_out_off}
  start_array['5SS']=${cov_in_exon}
  unset end_array; declare -A end_array
  end_array['3SS']=${cov_in_exon}
  end_array['5SS']=${cov_out_off}

  # check for computation of ASE and CE lists
  local prefixes="${exon_table_prefixes_ref},${exon_table_prefixes}"
  local RRBS_filtering_prefixes=${exon_table_prefixes}
  if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    local RRBS_filtering_prefixes="${prefixes}"
  fi



  for SS in '3SS' '5SS'; do
    for RRBS_window in $( echo ${RRBS_windows} | tr ',' ' ' ); do
      # filter CN files for splicing regions and computing distance to splicing site (stranded), and compute the fractions of methCyt
      local cn_distToSS_dir=${res_dir}/CN_filtered/inEx${cov_in_exon}_outOff${cov_out_off}/distToSS/win${RRBS_window}b; mkdir -vp ${cn_distToSS_dir}
      local fract_out_dir=${cn_distToSS_dir}/fract_$( echo ${mC_rate_limits} | tr ',' '_' ); mkdir -vp ${fract_out_dir}
      echo ">>> filtering, distToSS and fraction" >&2
      local ss_bed_dir=${ss_dir}/win${RRBS_window}b/bed6
      cmd="
      RRBS_CNfiltering ${RRBS_filtering_prefixes} ${SS} ${exons_bed_dir} ${ss_bed_dir} ${mC_rate_limits} ${RRBS_cn_files} ${cn_distToSS_dir} ${fract_out_dir} ${start_array[$SS]} ${end_array[$SS]} ${RRBS_window} ${RRBS_cond} ${RRBS_rep} \
        --nbr-cores ${NBRCORES} ${with_refs_exons_flag}
        "
      # echo "$cmd"
      eval "$cmd"

      # build the figures
      echo ">>> metagene" >&2
      unset comp_pair_array; declare -A comp_pair_array
      comp_pair_array['raw']=$( echo $comp_pair | cut -d ',' -f 1 )
      comp_pair_array['diff']=$comp_pair
      for mode in 'raw' 'diff'; do
        echo ">>> $mode" >&2
        local out_fig_dir=${fig_dir}/RRBS/inEx${cov_in_exon}_outOff${cov_out_off}/${mode}MethCyt/win${RRBS_window}b
        local out_fig_prefix=${out_fig_dir}/${fig_prefix}_${SS}
        cmd="
        RRBS_metagene ${prefixes} ${SS} ${out_fig_prefix} ${fract_out_dir} ${comp_pair_array[$mode]} ${mode} \
          ${color_pallette}
          "
        # echo "$cmd"
        eval "$cmd"
      done
    done
  done
}
