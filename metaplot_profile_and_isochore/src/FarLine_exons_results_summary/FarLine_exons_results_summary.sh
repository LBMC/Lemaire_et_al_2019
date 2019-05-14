#!/bin/bash

## Set the custom data
project_dir=$1 #~ /home/sebastien/work_temp/id_card_temp
CONF_FILE=$2 #~ ${project_dir}/configuration.sh
FarLine_analysis_dir=$3 #~ /home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/RNA-Seq_results/analyse_stat/analyse_stat_siDNMT3b-siGL2_median_paired_psi_corrigee

# check for a defined figure prefix
if [[ "x${@}" == x*"--fig-prefix"* ]]; then
  exons_list_name="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
else
  exons_list_name="$( basename $FarLine_analysis_dir | sed s/analyse_stat_// | sed s/_median_paired_psi_corrigee// )"
fi

## Start the main part
# Set the general environment
FarLine_exons_results_summary_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${FarLine_exons_results_summary_dir}/src

## source ${src_dir}/array_value2index.sh
source ${src_dir}/env.sh ${project_dir} ${src_dir}
source ${src_dir}/nb_ref_lists_fc.sh

#####
## create simple skipped exon list table from FarLine all data table
if [[ "x${@}" == x*"--all-data-tab"* || "x${@}" == x*"--all-parts"* ]]; then
  # check for provided query all_data.tab (or equivalent with the necessary columns)
  if [[ "x${@}" == x*"--query-all-data-tab"* ]]; then
    query_all_data_tab=$( echo "$@" | awk -F '--query-all-data-tab ' '{ print $2 }' | cut -d ' ' -f 1 )
    query_all_data_tab_flag="--query-all-data-tab ${query_all_data_tab}"
    FarLine_analysis_dir='_'
  fi

  # check for already computed all_data_tab from FarLine analyses
  if [[ "x${@}" == x*"--FarLine-all-data-dir"* ]]; then
    FarLine_all_data_dir=$( echo "$@" | awk -F '--FarLine-all-data-dir ' '{ print $2 }' | cut -d ' ' -f 1 )
    FarLine_all_data_dir_flag="--FarLine-all-data-dir ${FarLine_all_data_dir}"
  fi


  echo -e "\n>>> Create the 'all data tables'." >&2
  source ${src_dir}/all_data_tab/all_data_tab_script.sh ${FarLine_analysis_dir} ${res_dir} ${exons_list_name} ${GENE_FDB_TAB} ${EXON_FDB_TAB} \
    ${query_all_data_tab_flag} ${FarLine_all_data_dir_flag}
fi
# exit

#####
## build the figure of results for each table of exons
if [[ "x${@}" == x*"--skipped-exon-list-results-summary"* || "x${@}" == x*"--all-parts"* ]]; then
  # check for results of reference exons
  # with_refs_exons_flag="--refs-exons-dir ${results_ASE_CE_dir}" #~ "--with-refs-exons"
  if [[ "x${@}" == x*"--no-ref-links"* ]]; then
    no_ref_links="--no-ref-links"
  elif [[ "x${@}" == x*"--refs-exons-dir"* ]]; then
    with_refs_exons_flag="--refs-exons-dir $( echo "$@" | awk -F '--refs-exons-dir ' '{ print $2 }' | cut -d ' ' -f 1 )"
  elif [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    with_refs_exons_flag="--with-refs-exons"
  else
    echo "Specify if you want to compute results for reference exons: '--with-refs-exons', or the path to the precomputed results for reference exons: --refs-exons-dir <dir>" >&2
    exit
  fi

  # check in custom reference exons are specified (ASE/CE from FasterDB are the default)
  exon_tables_ref=${data_dir}/ref_exon_lists/ASE_FDB_exon_gene_pos.txt,${data_dir}/ref_exon_lists/CE_FDB_exon_gene_pos.txt
  exon_table_prefixes_ref='ASE_FDB','CE_FDB'
  if [[ "x${@}" == x*"--cust-ref-lists"* ]]; then
    exon_tables_ref="$( echo "$@" | awk -F '--cust-ref-lists ' '{ print $2 }' | cut -d ' ' -f 1 )"
    exon_table_prefixes_ref=$( echo "$@" | awk -F '--cust-ref-lists ' '{ print $2 }' | cut -d ' ' -f 2 )
    exon_tables_ref_array=( $( echo "$exon_tables_ref" | tr ',' ' ' ) )
    exon_tables_ref="--ex-tab-ref ${exon_tables_ref}"
  fi

  # check for specific colors in ggplots
  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    color_pallette="$( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
    color_pallette_array=( $( echo "$color_pallette" | tr ',' ' ' ) )
  fi

  # set which part to run
  selrs_parts="--all"
  if [[ "x${@}" == x*"--selrs-parts"* ]]; then
    selrs_parts_flag="--selrs-parts $( echo ${@} | awk -F '--selrs-parts ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # check if coverage summary metaplot should be "scale harmonised"
  sh_arg=''
  if [[ "x${@}" == x*"--scaleHarm"* ]]; then
    sh_arg="--scaleHarm"
  fi

  # recover the sub-selection of coverage tracks
  coverage_summary_bank_sel_arg=''
  flag='--cov-bank-sel'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # get subparts of base fraction if specified
  bf_arg='--bf-parts --bf-all'
  flag='--bf-parts'
  if [[ "x${@}" == x*"${flag}"* ]]; then
    bf_arg="--bf-parts $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # launch the computing
  echo ">>> SELRS for query exons" >&2

    # recover the path to the exon tables computed in all_data_tab part, and the corresponding color_pallette
  all_data_tab_dir=${res_dir}/all_data_tab
  exons_up_table=${all_data_tab_dir}/exon_up_${exons_list_name}.tab
  exons_down_table=${all_data_tab_dir}/exon_down_${exons_list_name}.tab

  unset color_pallette_ref_array; declare -a color_pallette_ref_array
  nbr=${#exon_tables_ref_array[@]}
  for ((xxx=0;xxx<${nbr};xxx++)); do
    color_pallette_ref_array=( ${color_pallette_ref_array[@]} ${color_pallette_array[$xxx]} )
  done

  echo -e "\n>>> Build the plots based on one table of skipping exons/events" >&2
  for exons_table in ${exons_up_table} ${exons_down_table}; do
    nana="$( tail -n +2 ${exons_table} | head -1 )"
    if [[ -n "$nana" ]]; then
      exons_tables="${exons_tables} ${exons_table}"
      color_pallette_ref_array=( ${color_pallette_ref_array[@]} ${color_pallette_array[$nbr]} )
    fi
    ((nbr++))
  done
  exons_tables=$( echo ${exons_tables} | tr ' ' '+' )
  color_pallette="--color-pallette $( echo ${color_pallette_ref_array[*]} | tr ' ' ',' )"

  # recover the path to the FarLine all data for psi-box ( with all the evaluated exons )
  # FarLine_all_data_psi_box=/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_methChromRegulators-siPP-siGL2/RNA-Seq_results/analyse_stat/analyse_stat_siPP-siGL2_median_paired_psi_corrigee/exon_skipping/tmp/all_data.tab
  FarLine_all_data_psi_box="--FL-adt-psibox $( echo ${@} | awk -F '--FarLine-all-data-psi-box ' '{ print $2 }' | cut -d ' ' -f 1 )"


  if [[ "x${@}" == x*"--no-main-selrs"* ]]; then
    echo ">>> Skip main SELRS computing on query exons." >&2
  else
    cmd="
    bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
      ${with_refs_exons_flag} ${FarLine_all_data_psi_box} --fig-prefix ${exons_list_name} ${color_pallette} ${selrs_parts_flag} ${no_ref_links} ${coverage_summary_bank_sel_arg} ${sh_arg} ${bf_arg}
    "
    echo "$cmd"
    eval "$cmd"
  fi

  ## check for analysis of neighbour exons
  for neighbour in 'up' 'down'; do
    if [[ "x${@}" == x*"--nb-${neighbour}"* ]]; then
      echo ">>> SELRS for ${neighbour}stream exons" >&2
      marque="nb${neighbour}"

      ## look for directory of precomputed results
      if [[ "x${@}" == x*"--refs-exons-dir"* ]]; then
        with_refs_exons_flag="--refs-exons-dir $( echo "$@" | awk -F "--nb-${neighbour} " '{ print $2 }' | cut -d ' ' -f 1 )"
      fi

      ## create the result dir for neighbour exons
      res_nb_dir="${res_dir}_${marque}"; mkdir -vp ${res_nb_dir}

      ## compute lists of neighbours of the ref exons
      exon_tables_ref_nb="${res_nb_dir}/$( echo "${exon_table_prefixes_ref}" | sed s%','%"_${marque}.tsv,${res_nb_dir}/"%g )_${marque}"
      exon_table_prefixes_ref_nb="$( echo "${exon_table_prefixes_ref}" | sed s/','/"_${marque},"/g )_${marque}"
      if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
        nb_ref_lists_fc ${exon_tables_ref} ${exon_tables_ref_nb} ${marque}
      fi

      ## compute the all_data_tabs for neighbours of the query exons ( 'exon_up' and 'exon_down' )
      exons_list_name_nb="${exons_list_name}_${marque}"
      exons_tables_nb="$( echo "${exons_tables}" | sed s,"${res_dir}","${res_nb_dir}",g | sed s/'\.tab'/"_${marque}.tab"/g )"
      Rscript ${src_dir}/nb_query_lists.r ${exons_tables} ${exons_tables_nb} ${marque} ${EXON_FDB_TAB}

      ## launch selrs with the replaced arguments
      cmd="
      bash ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh ${exons_tables_nb} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_nb_dir} ${FarLine_all_data_psi_box} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref_nb} ${exon_table_prefixes_ref_nb} \
        ${with_refs_exons_flag} ${no_ref_links} --fig-prefix ${exons_list_name_nb} ${color_pallette} ${selrs_parts_flag} ${sh_arg}
      "
      echo "$cmd"
      eval "$cmd"
    fi
  done

fi


#####
## build the resuming powerpoint
if [[ "x${@}" == x*"--pptx"* || "x${@}" == x*"--all-parts"* ]]; then
  bash ${src_dir}/pptx_builder/pptx_builder.sh ${res_dir}/figures ${project_dir}/${exons_list_name}_summary ${exons_list_name} ${res_dir} ${CONF_FILE}
fi
