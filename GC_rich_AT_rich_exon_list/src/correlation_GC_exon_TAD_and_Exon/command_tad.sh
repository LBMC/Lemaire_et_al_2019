#!/bin/bash

python3 src/correlation_GC_exon_TAD_and_Exon/create_GC_AT_bed_exon.py
grep -v "'MFE_3SS': None" result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons.bed > result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed
python3 src/correlation_GC_exon_TAD_and_Exon/tad_exon.py

# regulated exons bed
python3 src/correlation_GC_exon_TAD_and_Exon/branch_point_homogeneity_in_tad.py -t data/K562_Lieberman-raw_TADs.hg19.nochr.bed -e result/correlation_GC-AT-exons_TAD/data_for_regulated_exons.bed -o result/correlation_GC-AT-exons_TAD/ -n table_TAD_K562_regulated_exons


# k562
# CCE regultaed exons
python3 src/correlation_GC_exon_TAD_and_Exon/branch_point_homogeneity_in_tad.py -t data/K562_Lieberman-raw_TADs.hg19.nochr.bed -e result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed -o result/correlation_GC-AT-exons_TAD/ -n table_TAD_K562_regulated_CCE_exons


# MCF7
python3 src/correlation_GC_exon_TAD_and_Exon/branch_point_homogeneity_in_tad.py -t data/HiC_TADs_Stein-MCF7-WT.hg19.nochr.bed -e result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed -o result/correlation_GC-AT-exons_TAD/ -n table_TAD_MCF7_regulated_CCE_exons

