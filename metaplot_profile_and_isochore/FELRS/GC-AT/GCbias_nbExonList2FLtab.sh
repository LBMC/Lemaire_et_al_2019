#!/bin/bash

R_script=../../src/FarLine_exons_results_summary/src/misc_tools/exonList2FLtab_like.r

out_dir=./SFs_gathGCbiasCline_FLtab_nbh
mkdir -vp ${out_dir} >&2
echo ">>> ${out_dir}" >&2
for exon_list in $( find ./SFs_gathGCbiasCline_nbh/* -maxdepth 0 -name "*.txt" ); do
  echo "$exon_list" >&2
  exon_list_name=$( basename $exon_list .txt )
  cmd="Rscript ${R_script} ${exon_list} > ${out_dir}/${exon_list_name}.tab"
  # echo "$cmd"
  eval "$cmd"
done


out_dir=./SFs_gathGCbias_FLtab_nbh
mkdir -vp ${out_dir} >&2
echo ">>> ${out_dir}" >&2
for exon_list in $( find ./SFs_gathGCbias_nbh/* -maxdepth 0 -name "*.txt" ); do
  echo "$exon_list" >&2
  exon_list_name=$( basename $exon_list .txt )
  Rscript ${R_script} ${exon_list} > ${out_dir}/${exon_list_name}.tab
done


out_dir=./SFs_FLtab_nbh
mkdir -vp ${out_dir} >&2
echo ">>> ${out_dir}" >&2
for exon_list in $( find ./SFs_seb_nbh/* -maxdepth 0 -name "down*.txt" ); do
  echo "$exon_list" >&2
  exon_list_name=$( basename $exon_list .txt | sed s%down_regulated_exons_by_%%g | sed s%_in_%'_'%g | sed s%'.txt$'%%g )
  factor="$( echo $exon_list_name | sed s%'_.*$'%%g )"
  if [[ -n "$( grep "$factor" ~/work_temp/rna_folding_measures/data/GC_AT_rich_group_SFonly.tsv )" ]]; then
    cmd="Rscript ${R_script} ${exon_list} > ${out_dir}/${exon_list_name}.tab"
    # echo "$cmd"
    eval "$cmd"
  fi
done


out_dir=./ref_list_FLtab_nbh
mkdir -vp ${out_dir} >&2
echo ">>> ${out_dir}" >&2
for exon_list in $( find ./ref_list_nbh/* -maxdepth 0 -name "*.txt" ); do
  echo "$exon_list" >&2
  exon_list_name=$( basename $exon_list .txt )
  Rscript ${R_script} ${exon_list} > ${out_dir}/${exon_list_name}.tab
done


####
