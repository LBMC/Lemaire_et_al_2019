#!/bin/bash

python3 src/study_variance_exons_and_venn/variance_analysis.py
python3.5 src/study_variance_exons_and_venn/bed_file_creator_from_sed.py
mkdir result/variance_analysis/enrichment_barplots

python3.5 src/study_variance_exons_and_venn/check_enrichment_in_GC_AT_group.py -l result/variance_analysis/bed_file/small_sf-downregulated_exons_\(3-27nt\).bed result/variance_analysis/bed_file/small_CCE_exons_\(3-27nt\).bed result/variance_analysis/bed_file/CCE_exons_27nt+.bed -b Small_down_exons Small_cce_exons CCE_exons -o result/variance_analysis/enrichment_barplots
