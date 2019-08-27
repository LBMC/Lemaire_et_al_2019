#!/usr/bin/env bash

python3.5 src/study_microexons/Irimia/supplemental2bed.py

# Test
python3.5 src/study_microexons/Irimia/decorate_bed.py -i result/irimia_bed/test.bed -o result/irimia_bed/ -f data/fasterDB_lite.db -s data/sed.db -H data/Homo_sapiens.GRCh37.dna.primary_assembly.fa -v y

# Creating bed frequencies
python3.5 src/study_microexons/Irimia/decorate_bed.py -i result/irimia_bed/Irimia_et_al_microexons.bed -o result/irimia_bed/ -f data/fasterDB_lite.db -s data/sed.db -H data/Homo_sapiens.GRCh37.dna.primary_assembly.fa -v n
