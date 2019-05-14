#!/bin/bash

echo 'DADA' >&2

skipped_exon_list_results_summary_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

exon_table_list=$1
FEATURE_TAB=$2
FASTA_DIR=$3
CONF_FILE=$4
res_dir_selrss=$5
EXON_FDB_TAB=$6
CODING_EXON_TAB=$7
data_dir=$8
# exon_tables_ref=$9
# exon_table_prefixes_ref=${10}

# echo "$@" >&2

# check for computation of ASE and CE lists
if [[ "x${@}" == x*"--no-ref-links"* ]]; then
  echo "--- skip symbolic links to files of reference exons ---"
elif [[ "x${@}" == x*"--with-refs-exons"* ]]; then
  with_refs_exons_flag="--with-refs-exons"
else
  results_ASE_CE_dir="$( echo ${@} | awk -F '--refs-exons-dir ' '{ print $2 }' | cut -d ' ' -f 1 )"
  for yyy in $( echo ${results_ASE_CE_dir} | tr ',' ' ' );do
    for xxx in $( find "${yyy}/" -type f ); do
      lien_path=$( echo ${xxx} | sed s,"${yyy}","${res_dir_selrss}", )
      mkdir -vp $( dirname ${lien_path} )
      ln -fs ${xxx} ${lien_path}
    done
  done
  # connect to the existing ASE and CE results
fi

#change path to *_jbc directories
jbc_flag='';
if [[ "x${@}" == x*"--jbclaude"* ]]; then
    jbc_flag="--jbclaude"
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

# set which part to run
selrs_parts="--all"
if [[ "x${@}" == x*"--selrs-parts"* ]]; then
  selrs_parts="$( echo ${@} | awk -F '--selrs-parts ' '{ print $2 }' | cut -d ' ' -f 1 | tr ',' ' ' )"
fi

# get FarLine adt for deltaPSI 2D plot
flag='--FL-adt-psibox'
if [[ "x${@}" == x*"${flag}"* ]]; then
  FarLine_all_data_psi_box="--FL-adt-psibox $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 | tr ',' ' ' )"
fi

sh_arg=''
if [[ "x${@}" == x*"--scaleHarm"* ]]; then
  sh_arg="--scaleHarm"
fi

flag='--ex-tab-ref'
if [[ "x${@}" == x*"${flag}"* ]]; then
  exon_tables_ref="--ex-tab-ref $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  exon_table_prefixes_ref="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 2 )"
fi

coverage_summary_bank_sel_arg=''
flag='--cov-bank-sel'
if [[ "x${@}" == x*"${flag}"* ]]; then
  coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
fi

input_type_arg=''
if [[ "x${@}" == x*"--coord-only"* ]]; then
  input_type_arg="--coord-only"
fi

gene_concat_arg=''
if [[ "x${@}" == x*"--gene-concat"* ]]; then
  gene_concat_arg="--gene-concat"
fi

# get subparts of base fraction if specified
bf_arg='--bf-parts --bf-all'
flag='--bf-parts'
if [[ "x${@}" == x*"${flag}"* ]]; then
  bf_arg="--bf-parts $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
fi

various_length=''
flag='--various-length'
if [[ "x${@}" == x*"${flag}"* ]]; then
  various_length="--various-length"
fi

side_cut=''
flag='--side-cut'
if [[ "x${@}" == x*"${flag}"* ]]; then
  side_cut="--side-cut $( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
fi

echo "BOUH !" >&2

for exon_table in $( echo ${exon_table_list} | tr ',' ' ' ); do
  exon_table=$( echo $exon_table | tr '+' ',' )
  echo $exon_table >&2
  cmd="
  source ${skipped_exon_list_results_summary_script_dir}/skipped_exon_list_results_summary.sh ${exon_table} ${FEATURE_TAB} ${FASTA_DIR} ${CONF_FILE} ${res_dir_selrss} ${EXON_FDB_TAB} ${CODING_EXON_TAB} ${data_dir} ${exon_tables_ref} ${exon_table_prefixes_ref} \
    ${with_refs_exons_flag} ${FarLine_all_data_psi_box} --fig-prefix ${fig_prefix} ${color_pallette} ${coverage_summary_bank_sel_arg}\
    ${selrs_parts} ${sh_arg} ${input_type_arg} ${bf_arg} ${various_length} ${side_cut} ${gene_concat_arg} ${jbc_flag}\
    ;
    "
    # --all
    # --lists \
    # --features \
    # --psi-box \
    # --base-fraction \
    # --methCyt \
    # --coverage \

    echo "$cmd"
    eval "$cmd"
  # break
done
echo "FIN" >&2
exit

