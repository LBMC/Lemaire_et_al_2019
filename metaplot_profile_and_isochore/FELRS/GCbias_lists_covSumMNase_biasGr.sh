#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../src/FarLine_exons_results_summary/src

#### INPUT VARIABLES

anal_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/GC-AT'
all_res_dir=${anal_dir}/res
mkdir -vp "${all_res_dir}" >&2

exon_table_dir=${anal_dir}/SFs_gathGCbiasCline_FLtab #_gathGCbias_FLtab
ref_list_dir=${anal_dir}/ref_list

csvs_dir=${anal_dir}/../cov_sum_bank_sel
coverage_summary_bank_sel=all_MNase_tracks.txt

CONF_FILE='./configuration.sh'

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='--with-refs-exons' # should be no ref computing
no_ref_links="" #"--no-ref-links" # no ref in first place

color_pallette_GC="3361FF"
color_pallette_GC_wRd="66B2FF"
color_pallette_AT="006600"
color_pallette_AT_wRd="00BB00"

selrs_parts_flag='--selrs-parts --lists,--coverage' # pourrez rajouter '--features'
# selrs_parts_flag='--selrs-parts --coverage'

####

#### MAIN
## SELRS analysis
for nbh_sign in '' '_nbup' '_nbdown'; do
  if [[ "x$nbh_sign" == x"_nb"* ]]; then
    etd_suf='_nbh'
  else
    etd_suf=''
  fi

  # exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*" -type f | sort ) | tr ' ' '+' )
  exons_tables="${exon_table_dir}${etd_suf}/AT_rich${nbh_sign}.tab+${exon_table_dir}${etd_suf}/AT_rich_with_intersection${nbh_sign}.tab+${exon_table_dir}${etd_suf}/GC_rich${nbh_sign}.tab+${exon_table_dir}${etd_suf}/GC_rich_with_intersection${nbh_sign}.tab"
  echo "!!! $exons_tables" >&2
  color_pallette="--color-pallette ${color_pallette_AT},${color_pallette_AT_wRd},${color_pallette_GC},${color_pallette_GC_wRd}"

  ## ref lists
  ref_list_dir_run=${ref_list_dir}${etd_suf}

  ref_pref_array=( CCE ) #( all_intern ) #( CCE ACE )
  exon_table_prefixes_ref="$( echo "${ref_pref_array[*]}" | tr ' ' ',' )"
  for ref_pref_idx in ${!ref_pref_array[*]}; do
    ref_pref_array[$ref_pref_idx]="${ref_list_dir_run}/${ref_pref_array[$ref_pref_idx]}${nbh_sign}.txt"
  done
  exon_tables_ref="--ex-tab-ref $( echo "${ref_pref_array[*]}" | tr ' ' ',' )"

  color_pallette="--color-pallette CC0000,${color_pallette_AT},${color_pallette_AT_wRd},${color_pallette_GC},${color_pallette_GC_wRd}"


  # exons_list_name='BOUH'
  exons_list_name="GC-AT"
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  if [[ -n "$coverage_summary_bank_sel" ]]; then
    coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${coverage_summary_bank_sel} | tr ',' '\n' | awk -v csvs_dir=${csvs_dir} '{ ORS=","; print csvs_dir"/"$1 }' )"
  fi

  project_dir="${all_res_dir}/${exons_list_name}${nbh_sign}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg}
  "

  echo "$cmd" | tr ' ' '\n'
  eval "$cmd"
done


####
