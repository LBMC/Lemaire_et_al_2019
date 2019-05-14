#!/bin/bash

for input_dir in ref_list SFs_seb SFs_gathGCbias SFs_gathGCbiasCline; do
  # out_dir=SFs_gathGCbias_nbh
  # out_dir=SFs_seb_nbh
  out_dir=${input_dir}_nbh
  mkdir -vp $out_dir >&2

  dline=2
  # if [[ "$input_dir" == "SFs_gathGCbiasCline" || "$input_dir" == "SFs_gathGCbias" ]]; then
  #   dline=1
  # fi

  for exlist in $( ls ${input_dir}/*.txt ); do
    # exlist=SFs_gathGCbias/293T_AT_rich_group.txt

    out_file=${out_dir}/$( basename $exlist .txt )_nbup.txt
    if [[ "$dline" -eq "2" ]]; then head -1 $exlist > ${out_file}; fi
    tail -n +${dline} $exlist | awk -F '_' '{ OFS="_"; print $1,$2-1 }' >> $out_file

    out_file=${out_dir}/$( basename $exlist .txt )_nbdown.txt
    if [[ "$dline" -eq "2" ]]; then head -1 $exlist > ${out_file}; fi
    tail -n +${dline} $exlist | awk -F '_' '{ OFS="_"; print $1,$2+1 }' >> $out_file

  done
done

####
