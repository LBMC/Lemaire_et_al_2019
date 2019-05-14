#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../src/FarLine_exons_results_summary/src

#### INPUT VARIABLES

anal_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/GC-AT'
all_res_dir=${anal_dir}/res
mkdir -vp "${all_res_dir}" >&2

# group_tab=~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv
exon_table_dir=${anal_dir}/ref_list_FLtab #_gathGCbias_FLtab
# ref_list_dir=${anal_dir}/ref_list

csvs_dir=${anal_dir}/../cov_sum_bank_sel
coverage_summary_bank_sel=all_MNase_tracks.txt,good_ChIP_H3.txt

CONF_FILE='./configuration.sh'

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='--with-refs-exons' # should be no ref computing
no_ref_links="" #"--no-ref-links" # no ref in first place


selrs_parts_flag='--selrs-parts --lists,--base-fraction,--coverage' # pourrez rajouter '--features'
selrs_parts_flag='--selrs-parts --base-fraction,--coverage' # pourrez rajouter '--features'
bf_parts_flag='--bf-parts --dna-seq,--seq-to-fasta,--base-metaplot'


####

#### MAIN
## SELRS analysis
for nbh_sign in ''; do # '_nbup' '_nbdown'
  if [[ "x$nbh_sign" == x"_nb"* ]]; then
    etd_suf='_nbh'
  else
    etd_suf=''
  fi

  # exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*" -type f | sort ) | tr ' ' '+' )
  # exons_tables="${exon_table_dir}${etd_suf}/AT_rich_group${nbh_sign}.tab+${exon_table_dir}${etd_suf}/GC_rich_group${nbh_sign}.tab"
  exons_tables="${exon_table_dir}${etd_suf}/CCE${nbh_sign}.tab"

  # color_pallette="--color-pallette FF8000,CC0000,${color_pallette_AT},${color_pallette_GC}"
  color_pallette="--color-pallette CC0000"


  exons_list_name="CCE"
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${coverage_summary_bank_sel} | tr ',' '\n' | awk -v csvs_dir=${csvs_dir} '{ ORS=","; print csvs_dir"/"$1 }' )"

  project_dir="${all_res_dir}/${exons_list_name}${nbh_sign}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  cmd="
  bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
  ${with_refs_exons_flag} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg} ${bf_parts_flag}
  "

  echo "$cmd" | tr ' ' '\n'
  eval "$cmd"

  break
done


####
