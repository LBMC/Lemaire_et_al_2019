#!/bin/bash

expSeq_bw_dir=${data_dir}/bw_files/MBD-Seq
expSeq_bw_files="$( echo $( find ${expSeq_bw_dir}/ -maxdepth 1 -name "MBDseq_pool1_MCF7siGL2n?_*_treat_pileup.bw" | sort ) $( find ${expSeq_bw_dir}/ -maxdepth 1 -name "MBDseq_pool1_MCF7siPPn?_*_treat_pileup.bw" | sort ) | tr ' ' ',' )"
expSeq_cond=$( echo $( echo ${expSeq_bw_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 1 | sed s/'MBDseq_pool1_MCF7'/''/ ) | tr ' ' ',' )
expSeq_rep=$( echo $( echo ${expSeq_bw_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 2 | cut -d '_' -f 1 ) | tr ' ' ',' )
comp_pair=$( echo ${expSeq_cond} | tr ',' '\n' | uniq | tr '\n' ',' | sed s/','$/''/ )
hm_ymax=100 #maximum value in heatmap of coverage
mean_ymax=100
median_ymax=100
