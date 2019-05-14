#!/bin/bash

expSeq_bam_dir=${data_dir}/bam_files/RRBS_MCF7
expSeq_bam_files="$( echo $( find ${expSeq_bam_dir}/ -maxdepth 1 -name "RRBS_pool1_MCF7siGL2n?.bam" | sort ) $( find ${expSeq_bam_dir}/ -maxdepth 1 -name "RRBS_pool1_MCF7siPPn?.bam" | sort ) | tr ' ' ',' )"
expSeq_cond=$( echo $( echo ${expSeq_bam_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 1 | sed s/'RRBS_pool1_MCF7'/''/ ) | tr ' ' ',' )
expSeq_rep=$( echo $( echo ${expSeq_bam_files} | tr ',' '\n' | awk -F '/' '{ print $NF }' | cut -d 'n' -f 2 | cut -d '.' -f 1 ) | tr ' ' ',' )
comp_pair=$( echo ${expSeq_cond} | tr ',' '\n' | uniq | tr '\n' ',' | sed s/','$/''/ )
