#!/bin/bash

exon_list_results_summary_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${exon_list_results_summary_dir}/base_fraction_metagene/base_fraction_metagene_script.sh
source ${exon_list_results_summary_dir}/coverage_summary/coverage_summary_script.sh
source ${exon_list_results_summary_dir}/RRBS/RRBS_script.sh

exon_table=$1
FEATURE_TAB=$2
FASTA_DIR=$3
CONF_FILE=$4
res_dir_selrs=$5
EXON_FDB_TAB=$6
CODING_EXON_TAB=$7
data_dir=$8
# exon_tables_ref=$9
# exon_table_prefixes_ref=${10}

source ${CONF_FILE}

exon_table_prefixes=$( echo $( echo ${exon_table} | tr ',' '\n' | awk -F '/' '{ print $NF }' | sed s/'.tab'$// ) | tr ' ' ',' )

fig_dir=${res_dir_selrs}/figures; mkdir -vp ${fig_dir}
exons_list_dir=${res_dir_selrs}/exons_lists; mkdir -vp ${exons_list_dir}
feat_tab_dir=${res_dir_selrs}/feature_tabs; mkdir -vp ${feat_tab_dir}
ss_list_dir=${res_dir_selrs}/splicing_site_regions

# check for computation of reference exon lists
if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
  with_refs_exons_flag="--with-refs-exons"
fi

# check for a defined figure prefix
fig_prefix='test'
if [[ "x${@}" == x*"--fig-prefix"* ]]; then
  fig_prefix="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
fi

# check for specific colors in ggplots
if [[ "x${@}" == x*"--color-pallette"* ]]; then
  color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
fi

# get FarLine adt for deltaPSI 2D plot
flag='--FL-adt-psibox'
if [[ "x${@}" == x*"${flag}"* ]]; then
  FarLine_all_data_psi_box="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 | tr ',' ' ' )"
fi

sh_arg=''
if [[ "x${@}" == x*"--scaleHarm"* ]]; then
  sh_arg="--scaleHarm"
fi

flag='--ex-tab-ref'
if [[ "x${@}" == x*"${flag}"* ]]; then
  exon_tables_ref="--ex-tab-ref $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  exon_table_prefixes_ref="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 2 )"
  arg_exon_table_prefixes_ref="--ex-tab-ref ${exon_table_prefixes_ref}"
fi

input_type_arg=''
if [[ "x${@}" == x*"--coord-only"* ]]; then
  input_type_arg="--coord-only"
fi

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

various_length=''
cov_bed_dir=${ss_list_dir}
flag='--various-length'
if [[ "x${@}" == x*"${flag}"* ]]; then
  various_length="--various-length"
  cov_bed_dir=${exons_list_dir}
fi

side_cut=''
flag='--side-cut'
if [[ "x${@}" == x*"${flag}"* ]]; then
  side_cut="--side-cut $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
fi

gene_concat_arg=''
if [[ "x${@}" == x*"--gene-concat"* ]]; then
  gene_concat_arg="--gene-concat"
fi



###############################
## create the lists of exons
if [[ "x$@" == x*"--all"* || "x$@" == x*"--lists"* ]]; then
  echo -e "\n>>> Create the lists and BED files of exons, and regions of splicing sites" >&2
  sl2e_selrs_parts_flag=''
  if [[ "x$@" == x*"--all"* || "x$@" == x*"--lists "* || "x$@" == x*"--lists-base-fraction"* ]]; then
    sl2e_selrs_parts_flag="${sl2e_selrs_parts_flag} --base-fraction"
  fi
  if [[ "x$@" == x*"--all"* || "x$@" == x*"--lists "* || "x$@" == x*"--lists-coverage "* ]]; then
    sl2e_selrs_parts_flag="${sl2e_selrs_parts_flag} --coverage"
  fi

  cmd="
  bash ${exon_list_results_summary_dir}/skipping_list2exon_list/skipping_list2exon_list.sh ${exons_list_dir} ${FEATURE_TAB} ${exon_table} ${EXON_FDB_TAB} ${CONF_FILE} ${CODING_EXON_TAB} ${feat_tab_dir} \
    ${with_refs_exons_flag} ${exon_tables_ref} ${exon_table_prefixes_ref} ${input_type_arg} ${sl2e_selrs_parts_flag}
  "
  echo "$cmd"
  eval "$cmd"
fi


###############################
if [[ -n "${exon_table_prefixes_ref}" ]]; then
  ref_feat_tab="--ex-tab-ref ${feat_tab_dir}/$( echo ${exon_table_prefixes_ref} | sed s?','?"_feat.tab,${feat_tab_dir}/"?g )_feat.tab"
fi

## build the feature boxplot for the table of skipped exons
if [[ "x$@" == x*"--all"* || "x$@" == x*"--features"* ]]; then
  if [[ "x${@}" == x*"--coord-only"* ]]; then
    echo "!!! Input type is in coordinate format, Features can not be computed ! Must be FL-like tab of exons"
  else
    echo -e "\n>>> Draw the exon feature plots." >&2
    fig_exon_feature_plots_dir=${fig_dir}/exon_feature; mkdir -vp ${fig_exon_feature_plots_dir}
    ef_names=$( echo $( echo ${exon_table} | tr ',' '\n' | awk -F '/' '{ print $NF }' | sed s/'.tab'$// ) | tr ' ' ',' )
    cmd="
    Rscript ${exon_list_results_summary_dir}/exon_feature_plots/exon_feature_plots.r ${FEATURE_TAB} ${exon_table} ${fig_exon_feature_plots_dir} ${CODING_EXON_TAB} \
      ${ref_feat_tab} ${exon_table_prefixes_ref} --name ${ef_names} --fig-prefix ${fig_prefix} ${color_pallette}
      "
    echo "$cmd"
    eval "$cmd"
  fi
fi


## build boxplot of PSI (siGL2, siPP)
if [[ "x$@" == x*"--all"* || "x$@" == x*"--psi-box"* ]]; then
  if [[ "x${@}" == x*"--coord-only"* ]]; then
    echo "!!! Input type is in coordinate format, PSI boxplots can not be computed ! Must be FL-like tab of exons"
  else
    echo -e "\n>>> Draw the PSI boxplots." >&2
    fig_FarLine_psi_boxplots_dir=${fig_dir}/FarLine_psi; mkdir -vp ${fig_FarLine_psi_boxplots_dir}
    fpb_names=$( echo $( echo ${exon_table} | tr ',' '\n' | awk -F '/' '{ print $NF }' | sed s/'.tab'$// ) | tr ' ' ',' )
    cmd="
    Rscript ${exon_list_results_summary_dir}/FarLine_psi_boxplots/FarLine_psi_boxplots.r ${exon_table} ${FarLine_all_data_psi_box} ${FEATURE_TAB} ${CODING_EXON_TAB} ${fig_FarLine_psi_boxplots_dir} \
      ${ref_feat_tab} ${exon_table_prefixes_ref} --name ${fpb_names} --fig-prefix ${fig_prefix} ${color_pallette}
      "
    echo "$cmd"
    eval "$cmd"
  fi
fi


## build the figure of base fraction metagenes
if [[ "x$@" == x*"--all"* || "x$@" == x*"--base-fraction"* ]]; then
  echo -e "\n>>> Draw the base fraction metagenes'.">&2
  cmd="
  base_fraction_metagene_script ${exon_table_prefixes} ${FASTA_DIR} ${res_dir_selrs} ${CONF_FILE} \
    ${with_refs_exons_flag} ${arg_exon_table_prefixes_ref} --fig-prefix ${fig_prefix} ${bf_arg}
    "
  echo "$cmd"
  eval "$cmd"
fi


# ## build metagenes of RRBS 5mC rate
# if [[ "x$@" == x*"--all"* || "x$@" == x*"--methCyt"* ]]; then
#   echo -e "\n>>> Draw the metagenes of methyl-cytosine rate"
#   cmd="
#   RRBS_script ${exon_table_prefixes} ${res_dir_selrs} ${data_dir} ${exons_list_dir}/bed6 ${ss_list_dir} ${CONF_FILE} ${fig_dir} ${exon_table_prefixes_ref} \
#     ${with_refs_exons_flag} --fig-prefix ${fig_prefix} ${color_pallette}
#     "
#   echo "$cmd"
#   eval "$cmd"
# fi


## build metagenes of sequencing coverages
if [[ "x$@" == x*"--all"* || "x$@" == x*"--coverage"* ]]; then
  echo -e "\n>>> Draw the heatmaps and metagenes of sequencing coverages'.">&2
  cmd="
  coverage_summary_script ${exon_table_prefixes} ${cov_bed_dir} ${fig_dir} ${CONF_FILE} ${data_dir} \
    ${with_refs_exons_flag} ${arg_exon_table_prefixes_ref} --fig-prefix ${fig_prefix} ${color_pallette} ${coverage_summary_bank_sel_arg} ${sh_arg} ${various_length} ${side_cut} ${gene_concat_arg}
    "
  echo "$cmd"
  eval "$cmd"
fi
