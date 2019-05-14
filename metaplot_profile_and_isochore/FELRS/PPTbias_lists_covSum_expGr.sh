#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../src/FarLine_exons_results_summary/src

#### INPUT VARIABLES

anal_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/GC-AT'
all_res_dir=${anal_dir}/res
mkdir -vp "${all_res_dir}" >&2

group_tab=~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv
color_tab=~/work_temp/rna_folding_measures/data/force_acceptor_union_dataset_wColor.txt
exon_table_dir=${anal_dir}/SFs_FLtab #_gathGCbias_FLtab

csvs_dir=${anal_dir}/cov_sum_bank_sel
unset coverage_summary_bank_sel_array; declare -gA coverage_summary_bank_sel_array
coverage_summary_bank_sel_array['Hela']=Hela_MNase_tracks.txt
coverage_summary_bank_sel_array['K562']=K562_MNase_tracks.txt

CONF_FILE='./configuration.sh'

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='' # should be no ref computing
no_ref_links="--no-ref-links" # no ref in first place

# color_pallette_GC="FF0000"
# color_pallette_AT="0088FF"

selrs_parts_flag='--selrs-parts --lists,--coverage' # pourrez rajouter '--features'
# selrs_parts_flag='--selrs-parts --coverage'

####

#### MAIN
## get list of splicing factor per cell line
sf_cline_list="$( ls ${exon_table_dir}/ | \
sed s%.tab%%g | \
sed -E s%_[-a-zA-Z0-9]*$%' &'%g | \
sed s%' _'%'\t'% | \
sort -k2,2 | uniq )"

cline_set="$( echo "$sf_cline_list" | \
cut -d "$( printf '\t' )" -f 2 | \
sort | uniq )"

unset sf_list; declare -gA sf_list
for cline in $cline_set; do
  sf_list[${cline}]="$( echo "$sf_cline_list" | \
  grep -E "${cline}$" | \
  awk '{ ORS=","; print $1 }' | \
  sed s%',$'%% )"
done


## SELRS analysis
## query exons

cline_list="$( echo ${!coverage_summary_bank_sel_array[*]} )"
# cline_list="Hela K562 MCF7"
for nbh_sign in '' '_nbup' '_nbdown'; do
  if [[ "x$nbh_sign" == x"_nb"* ]]; then
    etd_suf='_nbh'
  else
    etd_suf=''
  fi

  for cline in $( echo $cline_list | tr ',' ' ' ); do

    echo ">>> $cline" >&2

    # exons_tables=$( echo $( find ${exon_table_dir}/* -maxdepth 0 -name "${cline}_*" -type f | sort ) | tr ' ' '+' )
    color_pallette=''
    exons_tables=''
    for factor in $( echo ${sf_list[$cline]} | tr ',' ' ' ); do
      table1="$( find ${exon_table_dir}${etd_suf}/* -maxdepth 0 -name "${factor}_${cline}${nbh_sign}.tab" )"
      exons_tables="$exons_tables $table1"
      gc_group="$( grep $factor ${group_tab} | awk '{ print $2 }' )"
      # if [[ "x$gc_group" == x"AT"* ]]; then
      #   color="0088FF"
      # elif [[ "x$gc_group" == x"GC"* ]]; then
      #   color="FF0000"
      # fi
      color="$( grep -E "^$factor" ${color_tab} | awk '{ print $NF }' | sed s%'#'%''% )"
      color_pallette="${color_pallette} $color"
    done
    exons_tables="$( echo $exons_tables | tr ' ' '+' )"
    color_pallette="--color-pallette $( echo $color_pallette | tr ' ' ',' )"


    # exons_list_name='BOUH'
    exons_list_name="${cline}${nbh_sign}"
    echo ">>> EXON LIST NAME: $exons_list_name" >&2

    coverage_summary_bank_sel="${csvs_dir}/${coverage_summary_bank_sel_array[${cline}]}"
    coverage_summary_bank_sel_arg="--cov-bank-sel ${coverage_summary_bank_sel}"

    project_dir="${all_res_dir}/${exons_list_name}_PPT_bias_byExp"
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
done


####
