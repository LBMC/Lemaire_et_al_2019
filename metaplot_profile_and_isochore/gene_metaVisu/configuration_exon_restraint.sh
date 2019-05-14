#!/bin/bash

## general
NBRCORES=1 #1

## base fraction parameters
window_sizes=20,1 #20,4,1
in_exon_list=50 #,200
out_off_list=500 #,200

## mCpG metagene
min_reads=10

## coverage parameters
cov_in_exon=50 #200
cov_out_off=150 #200
cov_center_off=100
various_length_nbins=1000

## clustering parameters
nbr_clusters=5
min_max_signals=0
