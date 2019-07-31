#!/bin/bash

mkdir result/candidate
mkdir result/candidate/GC
python3.5 src/candidate_exons/GC_candidate_exons.py -e result/GC_rich_exons -n GC -l SNRNP70 SNRPC DDX5_DDX17 -s data/sed.db -f data/fasterDB_lite.db -o result/candidate/GC
python3.5 src/candidate_exons/GC_candidate_exons.py -e result/GC_rich_exons -n GC -l SNRPC DDX5_DDX17 -s data/sed.db -f data/fasterDB_lite.db -o result/candidate/GC

mkdir result/candidate/AT
python3.5 src/candidate_exons/GC_candidate_exons.py -e result/AT_rich_exons -n AT -l U2AF2 SF1 -s data/sed.db -f data/fasterDB_lite.db -o result/candidate/AT --ss "3'ss"

