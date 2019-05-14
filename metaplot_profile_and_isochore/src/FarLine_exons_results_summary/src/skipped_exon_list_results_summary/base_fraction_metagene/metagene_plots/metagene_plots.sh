#!/bin/bash

metagene_plots_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# source ${metagene_plots_dir}/base_fraction_metagene/start_rel.r

base_vec_arg_trans() {
  echo "$1" | sed s/'[ATGC]'/'&,'/g | sed s/',p'//g | sed s/','$//
}

metagene_plots () {
  # local CONF_FILE=$6
  local window_size=$5
  local in_exon=$6
  local out_off=$7

  local fig_dir=$4

  local exon_seq_files=$1
  local exon_seq_names=$2
  local SS_type=$3 #~ 3SS or 5SS

  local NBRCORES=1
  local nbr_proc=0
  local nbr_task=0
  if [[ "x$@" == x*"--nbr-cores"* ]]; then
    NBRCORES=$( echo "$@" | awk -F '--nbr-cores' '{ print $NF }' | cut -d ' ' -f 2 )
  fi

  local fig_prefix='test'
  if [[ "x$@" == x*"--fig-prefix"* ]]; then
    local fig_prefix=$( echo "$@" | awk -F '--fig-prefix' '{ print $NF }' | cut -d ' ' -f 2 )
  fi

  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi


  ## build the metàgene of %base from the table of sequences
  # set the relative positional start of the sequences
  # source $CONF_FILE
  local pos_index=1

  unset start_rel; declare -A start_rel
  local start_rel_script=${metagene_plots_dir}/base_fraction_metagene/start_rel.r
  local start_rel['3SS']=$( Rscript ${start_rel_script} "$(($pos_index-$out_off))" $window_size 1 )
  local start_rel['5SS']=$( Rscript ${start_rel_script} "-$(($pos_index+$in_exon-1))" $window_size 1 )
  # echo ${start_rel[*]} >&2

  unset ylim_array; declare -A ylim_array
  local ylim_array['A']='0.07,0.35'
  local ylim_array['T']='0.1,0.55'
  local ylim_array['G']='0.1,0.4'
  local ylim_array['C']='0.15,0.5'
  local ylim_array['AT']='0.25,0.75'
  local ylim_array['GC']='0.3,0.65'
  local ylim_array['AG']='0.2,0.65'
  local ylim_array['TC']='0.4,0.8'
  local ylim_array['AC']='0.3,0.7'
  local ylim_array['TG']='0.3,0.7'
  local ylim_array['CpG']='0,0.15'
  local ylim_array['GpC']='0,0.25'
  local ylim_array['ApT']='0,0.25'
  local ylim_array['TpA']='0,0.25'
  local ylim_array['CpA']='0,0.25'
  local ylim_array['ApC']='0,0.15'
  local ylim_array['GpA']='0,0.25'
  local ylim_array['ApG']='0,0.25'
  local ylim_array['CpT']='0.1,0.35'
  local ylim_array['TpC']='0.05,0.3'
  local ylim_array['GpT']='0,0.25'
  local ylim_array['TpG']='0.1,0.25'
  local ylim_array['ApA']='0,0.2'
  local ylim_array['TpT']='0,0.4'
  local ylim_array['GpG']='0,0.25'
  local ylim_array['CpC']='0.05,0.35'
  local ylim_array['GpGpG']='0,0.25'
  local ylim_array['CpCpC']='0,0.25'
  local ylim_array['GGGAGG']='0,0.1'
  local ylim_array['G+C']='-0.7,0.4'
  local ylim_array['A+T']='-0.7,0.25'
  local ylim_array['C+T']='-0.7,0.3'
  local ylim_array['C+A']='-0.7,1'
  local ylim_array['semiG4-coding']='0,0.6'
  local ylim_array['truncatedG4-coding']='0,0.25'
  local ylim_array['G4-coding']='0,0.2'
  local ylim_array['G4-transcribed']='0,0.2'
  local ylim_array['A,T,G,C']='0,0.25'
  local ylim_array['AT,TG,GC,CA,AG,TC']='0,0.13'
  local ylim_array['CGT,AGT,ACT,ACG']='0,0.05'
  local ylim_array['GpTTpG']='0.15,0.4'
  local ylim_array['GpG']='0,0.4'
  local ylim_array['GpGpT']='0,0.25'
  local ylim_array['thinSemiG4-coding']='0,0.8'
  local ylim_array['thinG4-coding']='0,0.25'
  local ylim_array['relaxedSemiG4-coding']='0,0.75'
  local ylim_array['G6-1MM']='0,0.5'
  local ylim_array['RRACH']='0,0.15'
  local ylim_array['DGACH']='0,0.1'
  local ylim_array['DSAGR']='0,0.2'
  local ylim_array['GSARS']='0,0.12'
  local ylim_array['KSAGG']='0,0.15'
  local ylim_array['m6Asite']='0,0.3'

  # build the metagene for nucleotide nitrogen base content
  echo '---simple' >&2
  local base_array=( 'A' 'T' 'G' 'C' 'AT' 'GC' 'AG' 'TC' 'AC' 'TG' )
  for base_vec in ${base_array[*]}; do
    # build the metagenes
    # local base_vec_arg="$( echo "$base_vec" | sed s/'[ATGC]'/'&,'/ | sed s/',p'//g | sed s/','$// )"
    local base_vec_arg="$( base_vec_arg_trans "$base_vec" )"
    ylim_arg=''
    if [ -n "${ylim_array[${base_vec}]}" ]; then
      local ylim_arg="--ylim ${ylim_array[${base_vec}]}"
    fi
    Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
      "${exon_seq_files}" \
      "${exon_seq_names}" \
      $window_size \
      ${start_rel[$SS_type]} \
      ${fig_dir}/${fig_prefix} \
      "$base_vec_arg" \
      ${ylim_arg} \
      ${color_pallette} &

    nbr_proc=$(($nbr_proc+1))
    nbr_task=$(($nbr_task+1))
    if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
      nbr_proc=0
      wait
    fi

  done
  # wait

  # # build the metagene for pairs of nucleotide nitrogen base content
  # echo '---duo' >&2
  # local base_array=( 'CpG' 'GpC' 'ApT' 'TpA' 'CpA' 'ApC' 'GpA' 'ApG' 'CpT' 'TpC' 'GpT' 'TpG' 'ApA' 'TpT' 'GpG' 'CpC' )
  # for base_vec in ${base_array[*]}; do
  #   # build the metagenes
  #   # local base_vec_arg="$( echo "$base_vec" | sed s/'[ATGC]'/'&,'/ | sed s/',p'//g | sed s/','$// )"
  #   local base_vec_arg="$( base_vec_arg_trans "$base_vec" )"
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${base_vec}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #   nbr_proc=$(($nbr_proc+1))
  #   nbr_task=$(($nbr_task+1))
  #   if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #     nbr_proc=0
  #     wait
  #   fi
  #
  # done
  # # wait
  #
  # # build the metagene for G-quadruplexes properties
  # echo '---G or C triplets' >&2
  # local base_array=( 'GpGpG' 'CpCpC' )
  # for base_vec in ${base_array[*]}; do
  #   # build the metagenes
  #   # local base_vec_arg="$( echo "$base_vec" | sed s/'[ATGC]'/'&,'/g | sed s/',p'//g | sed s/','$// )"
  #   local base_vec_arg="$( base_vec_arg_trans "$base_vec" )"
  #   # echo "$base_vec_arg">&2
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${base_vec}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #   nbr_proc=$(($nbr_proc+1))
  #   nbr_task=$(($nbr_task+1))
  #   if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #     nbr_proc=0
  #     wait
  #   fi
  #
  # done
  # # wait
  #
  # # specific patterns
  # echo '---specific patterns' >&2
  # local base_array=( 'GGGAGG' )
  # local regexpr_name_array=( 'hnRNPH1' )
  # # for base_vec in ${base_array[*]}; do
  # for (( G4_index=0; G4_index<${#base_array[*]}; G4_index++ )); do
  #   # build the metagenes
  #   local base_vec_arg="${base_array[$G4_index]}"
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${base_vec}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     --regexpr ${regexpr_name_array[$G4_index]} \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #     nbr_proc=$(($nbr_proc+1))
  #     nbr_task=$(($nbr_task+1))
  #     if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #       nbr_proc=0
  #       wait
  #     fi
  #
  #   done
  #   # wait
  #
  # # build the metagene for nucleotide nitrogen base skewness
  # echo '---skew' >&2
  # local base_array=( 'G+C' 'C+A' 'A+T' 'C+T' )
  # for base_comb in ${base_array[*]}; do
  #   # separate the base of interest and the other basis to take in consideration
  #   local base_vec=$( echo $base_comb | cut -d '+' -f 1 )
  #   local counter_base_vec=$( echo $base_comb | cut -d '+' -f 2 )
  #
  #   # build the metagenes
  #   local base_vec_arg="$( echo "$base_vec" | sed s/'[ATGC]'/'&,'/ | sed s/','$// )"
  #   local counter_base_vec_arg="$( echo "$counter_base_vec" | sed s/'[ATGC]'/'&,'/ | sed s/','$// )"
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${base_comb}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     ${color_pallette} \
  #     ${ylim_arg} \
  #     --skewness ${counter_base_vec_arg} &
  #
  #   nbr_proc=$(($nbr_proc+1))
  #   nbr_task=$(($nbr_task+1))
  #   if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #     nbr_proc=0
  #     wait
  #   fi
  #
  # done
  # # wait
  #
  # # build the metagene for semi G-quadruplex content
  # echo '---semiG4-coding' >&2
  # local base_array=( 'GGG+([ATGC]{1,6}GGG+)+' 'GGG+([ATGC]{1,6}GGG+){1}([ATGC]{1,6}GGG+)+' 'GGG+([ATGC]{1,6}GGG+){2}([ATGC]{1,6}GGG+)+' 'CCC+([ATGC]{1,6}CCC+){2}([ATGC]{1,6}CCC+)+' )
  # local regexpr_name_array=( 'semiG4-coding' 'truncatedG4-coding' 'G4-coding' 'G4-transcribed' )
  # # for base_vec in ${base_array[*]}; do
  # for (( G4_index=0; G4_index<${#base_array[*]}; G4_index++ )); do
  #   # build the metagenes
  #   local base_vec_arg="${base_array[$G4_index]}"
  #   local regexpr_name=${regexpr_name_array[$G4_index]}
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${regexpr_name}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${regexpr_name}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     --regexpr ${regexpr_name} \
  #     --presence-only \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #     nbr_proc=$(($nbr_proc+1))
  #     nbr_task=$(($nbr_task+1))
  #     if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #       nbr_proc=0
  #       wait
  #     fi
  #
  #   done
  #   # wait
  #
  # # build the profile of entropy (based on ATGC bases)
  # echo '---entropy' >&2
  # local ent_base_array=( 'A,T,G,C' 'AT,TG,GC,CA,AG,TC' 'CGT,AGT,ACT,ACG' )
  # for ent_base_vec in ${ent_base_array[*]}; do
  #   ylim_arg=''
  #   echo $ent_base_vec >&2
  #   if [ -n "${ylim_array[${ent_base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${ent_base_vec}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #   "${exon_seq_files}" \
  #   "${exon_seq_names}" \
  #   $window_size \
  #   ${start_rel[$SS_type]} \
  #   ${fig_dir}/${fig_prefix} \
  #   "-" \
  #   ${color_pallette} \
  #   ${ylim_arg} \
  #   --entropy ${ent_base_vec} --struct-info
  # done
  #
  # # build the metagene for pairs of nucleotide nitrogen base content
  # echo '---G4-GGU-repeats' >&2
  # local base_array=( 'GpTTpG' 'GpG' 'GpGpT' )
  # for base_vec in ${base_array[*]}; do
  #   echo $base_vec >&2
  #   # build the metagenes
  #   local base_vec_arg="$( base_vec_arg_trans "$base_vec" )"
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${base_vec}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${base_vec}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #   nbr_proc=$(($nbr_proc+1))
  #   nbr_task=$(($nbr_task+1))
  #   if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #     nbr_proc=0
  #     wait
  #   fi
  #
  # done
  # # wait
  #
  # local base_array=( 'GG+([ATGC]{1,6}GG+)+' 'GG+([ATGC]{1,6}GGG+){2}([ATGC]{1,6}GGG+)+' '(GGG+([ATGC]{1,6}GG+)|GG+([ATGC]{1,6}GGG+))' '([ATGC]?[ATGC]GGGGG[ATGC]?|[ATGC]?G[ATGC]GGGG[ATGC]?|[ATGC]?GG[ATGC]GGG[ATGC]?|[ATGC]?GGG[ATGC]GG[ATGC]?|[ATGC]?GGGG[ATGC]G[ATGC]?|[ATGC]?GGGGG[ATGC][ATGC]?)' )
  # local regexpr_name_array=( 'thinSemiG4-coding' 'thinG4-coding' 'relaxedSemiG4-coding' 'G6-1MM' )
  # # for base_vec in ${base_array[*]}; do
  # for (( G4_index=0; G4_index<${#base_array[*]}; G4_index++ )); do
  #   # build the metagenes
  #   local base_vec_arg="${base_array[$G4_index]}"
  #   local regexpr_name=${regexpr_name_array[$G4_index]}
  #   echo $regexpr_name >&2
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${regexpr_name}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${regexpr_name}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     --regexpr ${regexpr_name_array[$G4_index]} \
  #     --presence-only \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #     nbr_proc=$(($nbr_proc+1))
  #     nbr_task=$(($nbr_task+1))
  #     if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #       nbr_proc=0
  #       wait
  #     fi
  #
  #   done
  #
  # # build the metagene for one consensus site of m6A (link to U-track bound by hnRNP-K)
  # echo '---m6A site pattern' >&2
  # local base_array=( '(G|A)(G|A)AC(A|C|T)' '(G|A|T)GAC(A|T|C)' '(G|A|T)(G|C)AG(G|A)' 'G(G|C)A(G|A)(G|C)' '(G|T)(G|C)AGG' '((G|A)(G|A)AC(A|C|T)|(G|A|T)GAC(A|T|C)|(G|A|T)(G|C)AG(G|A)|G(G|C)A(G|A)(G|C)|(G|T)(G|C)AGG)' )
  # local regexpr_name_array=( 'RRACH' 'DGACH' 'DSAGR' 'GSARS' 'KSAGG' 'm6Asite' )
  # # for base_vec in ${base_array[*]}; do
  # for (( G4_index=0; G4_index<${#base_array[*]}; G4_index++ )); do
  #   # build the metagenes
  #   local base_vec_arg="${base_array[$G4_index]}"
  #   local regexpr_name=${regexpr_name_array[$G4_index]}
  #   echo $regexpr_name >&2
  #   ylim_arg=''
  #   if [ -n "${ylim_array[${regexpr_name}]}" ]; then
  #     local ylim_arg="--ylim ${ylim_array[${regexpr_name}]}"
  #   fi
  #   Rscript ${metagene_plots_dir}/base_fraction_metagene/base_fraction_metagene.r \
  #     "${exon_seq_files}" \
  #     "${exon_seq_names}" \
  #     $window_size \
  #     ${start_rel[$SS_type]} \
  #     ${fig_dir}/${fig_prefix} \
  #     "$base_vec_arg" \
  #     --regexpr ${regexpr_name_array[$G4_index]} \
  #     ${ylim_arg} \
  #     ${color_pallette} &
  #
  #     nbr_proc=$(($nbr_proc+1))
  #     nbr_task=$(($nbr_task+1))
  #     if [[ "$nbr_proc" == "$NBRCORES" || "$nbr_task" == "${#base_array[@]}" ]]; then
  #       nbr_proc=0
  #       wait
  #     fi
  #
  #   done

    wait


}
