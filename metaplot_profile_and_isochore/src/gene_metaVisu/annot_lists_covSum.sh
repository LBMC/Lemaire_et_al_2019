#!/bin/bash


selrs_example_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
src_dir=${selrs_example_dir}/../FarLine_exons_results_summary/src

#### INPUT VARIABLES

# anal_dir="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu"
# all_res_dir=${anal_dir}/res
geneReg_dir=$7
all_res_dir=$8

# change results path to *_jbc
jbc_path='';
jbc_flag='';
if [[ "x${@}" == x*"--jbclaude"* ]]; then
    jbc_path="_jbc"
    jbc_flag="--jbclaude"
    all_res_dir=${all_res_dir}${jbc_path}
elif [[ "x${@}" == x*"--alapendry"* ]]; then
    jbc_path="_lap"
    jbc_flag="--alapendry"
    all_res_dir=${all_res_dir}${jbc_path}
fi

if [[ "x${@}" == x*"--hg38"* ]]; then
  all_res_dir=${all_res_dir}_hg38
elif [[ "x${@}" == x*"--hg18"* ]]; then
  all_res_dir=${all_res_dir}_hg18
fi
mkdir -vp "${all_res_dir}" >&2

# group_tab=~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv
# exon_table_dir=${anal_dir}/SFs_gathGCbiasCline_FLtab #_gathGCbias_FLtab
# reg_table_dir=${anal_dir}/geneReg

# csvs_dir=${anal_dir}/csvs
# coverage_summary_bank_sel=$3 #POLR2A_tracks.txt #ARG
# coverage_summary_bank_sel_arg="--cov-bank-sel $( echo ${coverage_summary_bank_sel} | sed s%'[0-9a-zA-Z\._-]*'%"${csvs_dir}/&"%g)"
coverage_summary_bank_sel_arg="--cov-bank-sel $3"

exon_tables_ref='' # no ref (or PSIref with expressed exons)
exon_table_prefixes_ref='' # empty

with_refs_exons_flag='' # should be no ref computing
no_ref_links="--no-ref-links" # no ref in first place

color_pallette="--color-pallette $2" #ARG
# color_pallette_GC="FF0000"
# color_pallette_AT="0088FF"

selrs_parts_flag="--selrs-parts $6" #'--selrs-parts --lists,--coverage'
# selrs_parts_flag='--selrs-parts --lists'
# selrs_parts_flag='--selrs-parts --coverage'

####

#### MAIN
## SELRS analysis
annot_type_list=( $( echo $4 | tr ',' ' ' ) ) #ARG
for annot_type in ${annot_type_list[*]} ; do
  echo "$annot_type" >&2
  if [[ "x$annot_type" == x"gene" || "x$annot_type" == x"exon_sync" || "x$annot_type" == x"exon_n"[mp]*"_sync" || "x$annot_type" == x"exon_first" || "x$annot_type" == x"exon_intern" || "x$annot_type" == x"exon_intern_simple" || "x$annot_type" == x"exon_last" || "x$annot_type" == x"intron" || "x$annot_type" == x"intron_n"[mp]*"_sync" ]]; then
    plot_type='--various-length'
  else
    plot_type=''
  fi

  extensions_arg=''
  CONF_FILE='./configuration_exon.sh'
  if [[ "x$annot_type" == x"gene" ]]; then
    extensions_arg='--ext-up 200 --ext-dw 200'
    CONF_FILE='./configuration_gene.sh'
  elif [[ "x$annot_type" == x"exon" || "x$annot_type" == x"exon_ssWinSmall" || "x$annot_type" == x"exon_n"[mp]* ]]; then
    CONF_FILE='./configuration_exon_restraint.sh'
  elif [[ "x$annot_type" == x"exon_ext" ]]; then
    CONF_FILE='./configuration_exon.sh'
  elif [[ "x$annot_type" == x"exon_lap"* ]]; then
    CONF_FILE='./configuration_exon_lap.sh'
  fi

  gene_concat=''
  side_cut=''
  if [[ "x$annot_type" == x"exon_intern" || "x$annot_type" == x"intron" ]]; then
    gene_concat='--gene-concat'
    side_cut='--side-cut 200,0'
  elif [[ "x$annot_type" == x"exon_intern_simple" ]]; then
    side_cut='--side-cut 200,0'
  fi

  exons_tables_dir=${geneReg_dir}/geneReg${jbc_path}_felrs
  if [[ "x${@}" == x*"--hg38"* ]]; then
    exons_tables_dir=${geneReg_dir}/geneReg${jbc_path}_hg38_felrs
  elif [[ "x${@}" == x*"--hg18"* ]]; then
    exons_tables_dir=${geneReg_dir}/geneReg${jbc_path}_hg18_felrs
  fi

  if [[ "x$annot_type" == x"exon_sync" || "x$annot_type" == x"exon_ssWinSmall" || "x$annot_type" == x"exon_lap" ]]; then
    annot_type="exon"
  elif [[ "x$annot_type" == x"exon_n"[mp]*"_sync" ]]; then
    annot_type="$( basename $annot_type _sync )"
  elif [[ "x$annot_type" == x"exon_intern_simple" ]]; then
    annot_type="exon_intern"
  elif [[ "x$annot_type" == x"intron_n"[mp]*"_sync" ]]; then
    annot_type="$( basename $annot_type _sync )"
  fi
  exons_tables="$( echo "$1" | sed s%'[0-9a-zA-Z_-]*'%"${exons_tables_dir}/&\_${annot_type}.tab"%g | tr ',' '+' )"
  # echo "!!! $exons_tables" >&2
  # color_pallette="--color-pallette ${color_pallette_AT},${color_pallette_GC}"


  # exons_list_name='BOUH'
  exons_list_name="${5}_${annot_type}" #"Nuc" #ARG
  echo ">>> EXON LIST NAME: $exons_list_name" >&2

  project_dir="${all_res_dir}/${exons_list_name}"
  echo ${project_dir}
  mkdir -vp ${project_dir} >&2

  source ${src_dir}/env.sh ${project_dir} ${src_dir}

  # à adapter pour metaplot_2pos (intern, intron, first, last) | normal (query, TSS, geneEnd)
  cmd="
  bash
  ${src_dir}/skipped_exon_list_results_summary/skipped_exon_list_results_summary_script.sh
  ${exons_tables}
  ${FEATURE_TAB}
  ${FASTA_DIR}
  ${CONF_FILE}
  ${res_dir}
  ${EXON_FDB_TAB}
  ${CODING_EXON_TAB}
  ${data_dir}
  ${exon_tables_ref}
  ${exon_table_prefixes_ref}
  ${with_refs_exons_flag}
  --fig-prefix ${exons_list_name}
  ${color_pallette}
  ${selrs_parts_flag}
  ${no_ref_links}
  ${coverage_summary_bank_sel_arg}
  --coord-only
  ${plot_type}
  ${side_cut}
  ${extensions_arg}
  ${gene_concat}
  ${jbc_flag}
  "

  echo "$cmd" | tr ' ' '\n'
  eval $cmd

done


####
