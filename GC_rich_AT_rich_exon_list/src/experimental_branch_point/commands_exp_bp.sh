#!/usr/bien/env bash

mkdir result/experimental_branch_point
python3 src/experimental_branch_point/branch_poin_bed_merged.py

python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
-b result/experimental_branch_point/mercer_2015_intron_end_filter.bed \
-n mercer_intron_filtered -a result/AT_rich_exons \
-g result/GC_rich_exons \
-f data/fasterDB_lite.db  \
-s data/sed.db \
-o result/experimental_branch_point


python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
-b result/experimental_branch_point/pineda_2018_5cov_filter_intron_end_filter.bed \
-n pineda_5cov_intron_filtered -a result/AT_rich_exons \
-g result/GC_rich_exons \
-f data/fasterDB_lite.db  \
-s data/sed.db \
-o result/experimental_branch_point


python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
-b result/experimental_branch_point/taggart_2017_5cov_filter_intron_end_filter.bed \
-n taggart_5cov_intron_filtered -a result/AT_rich_exons \
-g result/GC_rich_exons \
-f data/fasterDB_lite.db  \
-s data/sed.db \
-o result/experimental_branch_point


python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
-b result/experimental_branch_point/merged_branch_point_all_intron_filter.bed \
-n merged_intron_filtered -a result/AT_rich_exons \
-g result/GC_rich_exons \
-f data/fasterDB_lite.db  \
-s data/sed.db \
-o result/experimental_branch_point


python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
-b result/experimental_branch_point/pineda-taggart_filtered_merged.bed \
-n pineda-taggart_filtered -a result/AT_rich_exons \
-g result/GC_rich_exons \
-f data/fasterDB_lite.db  \
-s data/sed.db \
-o result/experimental_branch_point


python3 src/experimental_branch_point/create_predicted_bp_bed.py
python3 src/experimental_branch_point/compute_intersection_between_experiment_and_predicted_branch_point.py