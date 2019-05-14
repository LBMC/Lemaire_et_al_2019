#!/bin/bash

expSeq_bw_files="$( echo $( find ${expSeq_bw_dir}/ -maxdepth 1 -name "*.bw" | sort ) | tr ' ' ',' )"
expSeq_cond=$( echo ${expSeq_bw_files} | sed s/[^,]*/'public'/g )
expSeq_rep=$( echo ${expSeq_bw_files} | awk -F ',' '{ OFS=","; for(i=1;i<=NF;i++) { $i=i }; print }' )
comp_pair='public'
