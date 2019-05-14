#!/usr/bin/Rscript
require("bigWig", quietly=TRUE)
args <- commandArgs(trailingOnly = TRUE)
str(args)
print(args[1])
print(load.bigWig(args[1]))

#2018 11 JBC
# Tester si un bw utilise les refs chrX ou X
# Use as :
# bw= [my_bw_file]
# Rscript checkBigWig.R $bw > ${bw}.info
# test=$(grep -e chrom -A 1 ${bw}.info | grep -e chr[0-9XYM])
# if [ -z "$test" ] ; then
#    echo "ref as X"
#    chrflag="""
# else
#    echo "ref as chrX"
#    chrflag="--add-ref-chr"
# fi
#
