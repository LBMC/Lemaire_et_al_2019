#!/bin/bash

## Generate the PPTX reports of the Felrs analysis
script_dir="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/pdf"
measure_window=50,500 #50,200
center_window=1000

MNase_tracks_Didier='GSM971951_GSE39579_MNase_HEK293_B1s,GSE116770_HeLa_MNAse_B1s,GSE14022_HeLa_MNAse_B1s,GSE42951_HeLa_MNAse_B1s,GSE42951_HeLa_MNAse_B2s,GSE65644_HeLa_MNAse_B1s,GSE65644_HeLa_MNAse_B1s.t,GSM2083140_GSE78984_K562_MNAse_300U_B1s.t' #
H3_tracks_Didier_hg19='GSE97827_GSM2579061_HEK293_H3_B1s,GSE60526_HeLa_H3_B1s.t,GSM1846180_GSE71809_K562_H3_B1s.t' #
H3_tracks_Didier_hg38='GSM1482823_HEK293_H3_cistrome53314,GSM1918491_K562_H3_cistrome57381,GSM1918492_K562_H3_cistrome57382,GSM1918497_K562_H3_cistrome57387,GSM1918498_K562_H3_cistrome57388,GSM1369943_MCF7_H3_cistrome48593'

## specific variable setters
MNase_fc () {
  wdir="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/res"
  ref_list='geneVisu_exon'
  pptx_suffix=Nuc_${ref_list}
  measure_list=${MNase_tracks_Didier}
}

H3_hg19_fc () {
  wdir="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/res"
  ref_list='geneVisu_exon'
  pptx_suffix=Nuc_${ref_list}
  measure_list=${H3_tracks_Didier_hg19}
}

H3_hg38_fc () {
  wdir="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/res_hg38"
  ref_list='geneVisu_exon'
  pptx_suffix=Nuc_${ref_list}
  measure_list=${H3_tracks_Didier_hg38}
}


#### MAIN

echo ${script_dir}/pptx_measure_multiAnal_builder.py

MNase_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center_MNase \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_ss_MNase \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;

H3_hg19_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center_H3_hg19 \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_ss_H3_hg19 \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;

H3_hg38_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center_H3_hg38 \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_ss_H3_hg38 \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;

