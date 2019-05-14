#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../src/FarLine_exons_results_summary/src

####
anal_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/GC-AT'
all_res_dir=${anal_dir}/res
mkdir -vp "${all_res_dir}" >&2

group_tab=~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv

csvs_dir=${anal_dir}/../cov_sum_bank_sel
unset coverage_summary_bank_sel_array; declare -gA coverage_summary_bank_sel_array
coverage_summary_bank_sel_array['Hela']=Hela_MNase_tracks.txt
coverage_summary_bank_sel_array['K562']=K562_MNase_tracks.txt

####

#### MAIN
CONF_FILE='./configuration.sh'

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='' # should be no ref computing
no_ref_links="--no-ref-links" # no ref in first place

color_pallette_GC="0000FF"
color_pallette_AT="006600"
color_pallette="--color-pallette ${color_pallette_AT},${color_pallette_GC}" # GCrich: red, ATrich: blue

selrs_parts_flag='--selrs-parts --lists,--coverage' # pourrez rajouter '--features'
# selrs_parts_flag='--selrs-parts --coverage'


## SELRS analysis
## query exons
exon_table_dir=${anal_dir}/SFs_gathGCbias_FLtab

cline_list="$( echo ${!coverage_summary_bank_sel_array[*]} )"
# cline_list="Hela K562 MCF7"
for cline in $( echo $cline_list | tr ',' ' ' ); do

  echo ">>> $cline" >&2

  exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*" -type f | sort ) | tr ' ' '+' )
  # exons_tables='ATrich.tsv,GCrich.tsv' # GC-rich and AT-rich lists (in one cell line) (like in results/all_data_tab/exon_{up,down}_…)
  # exons_tables='./SELRS_example/bouh_1.tab+./SELRS_example/bouh_2.tab' # GC-rich and AT-rich lists (in one cell line) (like in results/all_data_tab/exon_{up,down}_…)

  # exons_list_name='BOUH'
  exons_list_name="$( echo "$exons_tables" | awk -F '_AT_rich_group.tab+' '{ print $1 }' | awk -F '/' '{ print $NF }' )"
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  coverage_summary_bank_sel="${csvs_dir}/${coverage_summary_bank_sel_array[${cline}]}"
  coverage_summary_bank_sel_arg="--cov-bank-sel ${coverage_summary_bank_sel}"

  project_dir="${all_res_dir}/${exons_list_name}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${FarLine_all_data_psi_box} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg}
  "

  echo "$cmd"
  eval "$cmd"

done


## upstream exon
exon_table_dir=${anal_dir}/SFs_gathGCbias_FLtab_nbh
for cline in $( echo $cline_list | tr ',' ' ' ); do
  echo ">>> $cline" >&2

  exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*_nbup.tab" -type f | sort ) | tr ' ' '+' )
  # exons_list_name='BOUH'

  exons_list_name="$( echo "$exons_tables" | awk -F '_AT_rich_group_nbup.tab+' '{ print $1 }' | awk -F '/' '{ print $NF }' )_nbup"

  coverage_summary_bank_sel="${csvs_dir}/${coverage_summary_bank_sel_array[${cline}]}"
  coverage_summary_bank_sel_arg="--cov-bank-sel ${coverage_summary_bank_sel}"

  project_dir="${all_res_dir}/${exons_list_name}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${FarLine_all_data_psi_box} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg}
  "

  echo "$cmd"
  eval "$cmd"

done


## downstream exon
exon_table_dir=${anal_dir}/SFs_gathGCbias_FLtab_nbh
for cline in $( echo $cline_list | tr ',' ' ' ); do
  echo ">>> $cline" >&2

  exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*_nbdown.tab" -type f | sort ) | tr ' ' '+' )

  exons_list_name="$( echo "$exons_tables" | awk -F '_AT_rich_group_nbdown.tab+' '{ print $1 }' | awk -F '/' '{ print $NF }' )_nbdown"
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  coverage_summary_bank_sel="${csvs_dir}/${coverage_summary_bank_sel_array[${cline}]}"
  coverage_summary_bank_sel_arg="--cov-bank-sel ${coverage_summary_bank_sel}"

  project_dir="${all_res_dir}/${exons_list_name}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${FarLine_all_data_psi_box} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg}
  "

  echo "$cmd"
  eval "$cmd"

done

####
