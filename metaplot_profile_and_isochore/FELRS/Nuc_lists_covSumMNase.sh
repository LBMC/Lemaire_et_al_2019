#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../src/FarLine_exons_results_summary/src

#### INPUT VARIABLES

anal_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/nucleosome'
all_res_dir=${anal_dir}/res
mkdir -vp "${all_res_dir}" >&2

# group_tab=~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv
# exon_table_dir=${anal_dir}/SFs_gathGCbiasCline_FLtab #_gathGCbias_FLtab
reg_table_dir=${anal_dir}/query

csvs_dir=${anal_dir}/cov_sum_bank_sel
# unset coverage_summary_bank_sel_array; declare -gA coverage_summary_bank_sel_array
# coverage_summary_bank_sel_array['Hela']=Hela_MNase_tracks.txt
# coverage_summary_bank_sel_array['K562']=K562_MNase_tracks.txt
coverage_summary_bank_sel=auboeuf_MCF7_MNase_tracks.txt

CONF_FILE='./configuration.sh'

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='' # should be no ref computing
no_ref_links="--no-ref-links" # no ref in first place

color_pallette="--color-pallette 000000"
# color_pallette_GC="FF0000"
# color_pallette_AT="0088FF"

selrs_parts_flag='--selrs-parts --lists,--coverage' # pourrez rajouter '--features'
selrs_parts_flag='--selrs-parts --lists' # pourrez rajouter '--features'
selrs_parts_flag='--selrs-parts --coverage' # pourrez rajouter '--features'
# selrs_parts_flag='--selrs-parts --coverage'

####

#### MAIN
## SELRS analysis
# for nbh_sign in '' '_nbup' '_nbdown'; do
  # if [[ "x$nbh_sign" == x"_nb"* ]]; then
  #   etd_suf='_nbh'
  # else
  #   etd_suf=''
  # fi

  # exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*" -type f | sort ) | tr ' ' '+' )
  exons_tables="${reg_table_dir}/siGL2_B1s_normPos_73-74.tab"
  echo "!!! $exons_tables" >&2
  # color_pallette="--color-pallette ${color_pallette_AT},${color_pallette_GC}"


  # exons_list_name='BOUH'
  exons_list_name="Nuc"
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  coverage_summary_bank_sel_arg="--cov-bank-sel ${csvs_dir}/${coverage_summary_bank_sel}"

  project_dir="${all_res_dir}/${exons_list_name}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg} --coord-only
  "

  echo "$cmd" | tr ' ' '\n'
  eval "$cmd"

# done


####
