#!/usr/bin/env bash

mkdir result/figure_2_3
mkdir result/figure_2_3/test
python3.5 src/figures_2_3_creator/figure_2_3_creator.py -l result/GC_rich_exons result/AT_rich_exons -n GC-exons AT-exons -s data/sed.db -f data/fasterDB_lite.db -o result/figure_2_3/test -S 'gc_rich_down' 'at_rich_down'


