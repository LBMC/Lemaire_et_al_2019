#!/bin/bash
# the commands in article_LEMAIRE_scripts/GC_rich_AT_rich_exon_list/src/study_microexons/Irimia/commands.sh must have been launched
#Those commands have beenlaunched from the location "article_LEMAIRE_scripts/GC_rich_AT_rich_exon_list"
python3 src/study_microexons/variance_analysis.py
python3.5 src/study_microexons/bed_file_creator_from_sed.py
mkdir result/variance_analysis/enrichment_barplots

python3.5 src/study_microexons/check_enrichment_in_GC_AT_group.py -l result/variance_analysis/bed_file/small_sf-downregulated_exons_\(3-27nt\).bed result/variance_analysis/bed_file/small_CCE_exons_\(3-27nt\).bed result/variance_analysis/bed_file/CCE_exons_27nt+.bed -b Small_down_exons Small_cce_exons CCE_exons -o result/variance_analysis/enrichment_barplots


# Micro-exons Irimia et al.
mkdir result/variance_analysis/enrichment_Irimia_barplots
python3.5 src/study_microexons/check_enrichment_in_GC_AT_group.py -l result/irimia_bed/Irimia_et_al_microexons_freq.bed result/variance_analysis/bed_file/CCE_exons_27nt+.bed -b Small_Irimia_exons CCE_exons -o result/variance_analysis/enrichment_Irimia_barplots
