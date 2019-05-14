#!/bin/bash

# metagene_coverage_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary'
metagene_coverage_dir='../src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary'
# bw_file_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data/bw_files/encode_Pol2'
bw_file_dir='../src/FarLine_exons_results_summary/data/bw_files/encode_Pol2'

out_dir_gnl='./GC-AT/res/GC-AT/'

for bw_file in $( ls ${bw_file_dir}/*/*.bw ); do
  # bw_file='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data/bw_files/encode_Pol2/POLR2A/ENCSR000BGD_ENCFF593AJA_POLR2A_2.bw'
  cond='public'
  rep='1'
  comp_pair='public'

  bed6_files='./AT_rich_group_gene.bed,./GC_rich_group_gene.bed'
  prefixes='AT-rich,GC-rich'

  sign='GC-AT_genes'
  out_dir=$( basename $bw_file .bw )/cov

  off_set=0

  color_pallette='FF0000'
  if [[ "$( basename $bw_file )" == 'ENC'* ]]; then
    ref_chr_arg='--ref-add-chr'
  else
    ref_chr_arg='' #'--ref-add-chr'
  fi

  various_length='--various-length 2000'
  extensions='--ext-up 200 --ext-dw 200'

  Rscript ${metagene_coverage_dir}/metagene_coverage.r \
  -bw ${bw_file} \
  -cond ${cond} \
  -rep ${rep} \
  -annot ${bed6_files} \
  -prefixes ${prefixes} \
  -sign ${sign} \
  -out_dir ${out_dir} \
  -off_set ${off_set} \
  -comp_pair ${comp_pair} \
  ${color_pallette} \
  ${hm_ymax} \
  ${mean_ymax} \
  ${median_ymax} \
  ${sh_arg} \
  ${ref_chr_arg} \
  ${various_length} \
  ${extensions} \
  ;
done


####
