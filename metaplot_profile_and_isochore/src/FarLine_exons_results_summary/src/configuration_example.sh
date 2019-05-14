#!/bin/bash

## general
NBRCORES=3

## base fraction parameters
window_sizes=20 #20,4,1
in_exon_list=50 #,200
out_off_list=500 #,200
center_off=1000

## mCpG metagene
min_reads=10

## coverage parameters
cov_in_exon=50 #200
cov_out_off=500 #200 #500
cov_center_off=1000

mC_rate_limits=13,87
RRBS_windows=20,10

## clustering parameters
nbr_clusters=5
min_max_signals=5
