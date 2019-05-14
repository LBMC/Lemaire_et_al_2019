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
cov_out_off=200 #200
cov_center_off=200
various_length_nbins=100
ext_up=50
ext_dw=50

## clustering parameters
nbr_clusters=5
min_max_signals=0
