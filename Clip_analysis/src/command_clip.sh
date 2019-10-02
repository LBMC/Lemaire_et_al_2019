#!/bin/bash 

mkdir data/input_oAT-GC/
sed "s/\t/_/g" ../GC_rich_AT_rich_exon_list/result/figure_2_3/exon_list/GC_unregulated_exons.txt > data/input_oAT-GC/GC_unregulated_exons.txt
sed "s/\t/_/g" ../GC_rich_AT_rich_exon_list/result/figure_2_3/exon_list/AT_unregulated_exons.txt > data/input_oAT-GC/AT_unregulated_exons.txt
mkdir data/input
python3 src/input_creator.py
python3 src/create_cce_file.py

mkdir result/clip_merged
mkdir result/clip_merged/coverage/
mkdir result/clip_merged/coverage_extended_GC_AT_exons/


python3 src/metaexon_coverage_ctrl.py -i data/merged_bed/ -f data/input/ -e data/exon.bed -c data/hg19.ren.chrom.sizes -o result/clip_merged/coverage/ -m ../metaplot_profile_and_isochore/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r --ctrl y && python3 src/coverage_figure_compacter.py -f result/clip_merged/coverage/figure/metagene_mean/ && montage -density 300 -geometry +1+1 -tile 1x8 -compress jpeg result/clip_merged/coverage/figure/metagene_mean/*_recap.pdf result/clip_merged/coverage/figure/metagene_mean/final_figures.pdf
python3 src/launcher_metaexon.py -i data/merged_bed/ -f data/input_oAT-GC/ -e data/exon.bed -c data/hg19.ren.chrom.sizes -o result/clip_merged/coverage_extended_GC_AT_exons/ -m ../metaplot_profile_and_isochore/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r && python3 src/coverage_figure_compacter.py -f result/clip_merged/coverage_extended_GC_AT_exons/figure/metagene_mean/ && montage -density 300 -geometry +1+1 -tile 1x8 -label %f -compress jpeg result/clip_merged/coverage_extended_GC_AT_exons/figure/metagene_mean/*_recap.pdf result/clip_merged/coverage_extended_GC_AT_exons/figure/metagene_mean/final_figures.pdf



