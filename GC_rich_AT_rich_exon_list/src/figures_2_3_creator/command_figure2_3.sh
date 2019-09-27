#!/usr/bin/env bash

mkdir result/figure_2_3
python3.5 src/figures_2_3_creator/custom_exon_list_creator.py
mkdir result/figure_2_3/test
python3.5 src/figures_2_3_creator/figure_2_3_creator.py -l result/GC_rich_exons result/AT_rich_exons -n GC-exons AT-exons -s data/sed.db -f data/fasterDB_lite.db -o result/figure_2_3/test2 -b ../Clip_analysis/data/coverage_project_selected/ -r ../Clip_analysis/data/hg19.ren.chrom.sizes -m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r -S 'gc_rich_down' 'at_rich_down'


mkdir result/figure_2_3/GA_CT_figures
mkdir result/figure_2_3/otherGC_GC_figures
mkdir result/figure_2_3/otherGC_AT_figures
mkdir result/figure_2_3/unregulated_GC_AT_exons
mkdir result/figure_2_3/otherGC_AT_GC_figures

python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
-l result/figure_2_3/exon_list/CT_rich_exons.txt result/figure_2_3/exon_list/GA_rich_exons.txt \
-n CT-exons GA-exons \
-N R \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/GA_CT_figures \
-b ../Clip_analysis/data/coverage_project_selected/ \
-r ../Clip_analysis/data/hg19.ren.chrom.sizes \
-m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r \
-S 'ct_rich_down' 'ga_rich_down' \
--reverse y


python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
-l result/figure_2_3/exon_list/GC_rich_exons.txt result/figure_2_3/exon_list/other_GC_rich_exons.txt \
-n GC-exons-bis otherGC-exons \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/otherGC_GC_figures \
-b ../Clip_analysis/data/coverage_project_selected/ \
-r ../Clip_analysis/data/hg19.ren.chrom.sizes \
-m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r \
-S 'gc_rich_down' 'other'


python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
-l result/figure_2_3/exon_list/other_GC_rich4ATcomp_exons.txt result/figure_2_3/exon_list/AT_rich_exons.txt \
-n otherGC-exons2 AT-exons-bis \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/otherGC_AT_figures \
-b ../Clip_analysis/data/coverage_project_selected/ \
-r ../Clip_analysis/data/hg19.ren.chrom.sizes \
-m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r \
-S 'other' 'at_rich_down'


python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
-l result/GC_rich_exons result/AT_rich_exons result/figure_2_3/exon_list/other_GC_exons/SRSF2_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/SRSF3_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/HNRNPC_all_exons.txt \
-n GC-exons AT-exons SRSF2-down-all SRSF3-down-all hnRNPC-down-all \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/otherGC_AT_GC_figures \
-b ../Clip_analysis/data/coverage_project_selected/ \
-r ../Clip_analysis/data/hg19.ren.chrom.sizes \
-m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r \


python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
-l result/figure_2_3/exon_list/GC_unregulated_exons.txt result/figure_2_3/exon_list/AT_unregulated_exons.txt \
-n GC-unreg-exons AT-unreg-exons \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/unregulated_GC_AT_exons \
-b ../Clip_analysis/data/coverage_project_selected/ \
-r ../Clip_analysis/data/hg19.ren.chrom.sizes \
-m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r


# Creation of custom figures
mkdir result/figure_2_3/custom_figure_oGC_exons
python3.5 src/figures_2_3_creator/custom_figure_creator.py \
-l result/GC_rich_exons result/AT_rich_exons result/figure_2_3/exon_list/other_GC_exons/SRSF2_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/SRSF3_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/HNRNPC_all_exons.txt \
-n GC-exons AT-exons SRSF2-down-all SRSF3-down-all hnRNPC-down-all \
-s data/sed.db \
-f data/fasterDB_lite.db \
-o result/figure_2_3/custom_figure_oGC_exons \
-t force_donor force_acceptor iupac_exon upstream_intron_size downstream_intron_size exon_size
