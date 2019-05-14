#!/bin/bash

expSeq_bw_dir=${data_dir}/bw_files/MNase-Seq
expSeq_bw_files="$( echo $( find ${expSeq_bw_dir}/ -maxdepth 1 -name "siGL2n1_treat_pileup.bw" | sort ) $( find ${expSeq_bw_dir}/ -maxdepth 1 -name "siDDX5-17n1_treat_pileup.bw" | sort ) | tr ' ' ',' )"
expSeq_cond=$( echo $( echo ${expSeq_bw_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 1 ) | tr ' ' ',' )
expSeq_rep=$( echo $( echo ${expSeq_bw_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 2 | cut -d '_' -f 1 ) | tr ' ' ',' )
comp_pair=$( echo ${expSeq_cond} | tr ',' '\n' | uniq | tr '\n' ',' | sed s/','$/''/ )
hm_ymax=100
mean_ymax=73
median_ymax=100
