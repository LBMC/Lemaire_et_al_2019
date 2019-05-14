#!/bin/bash

input_res_dir=$1 # res/geneVisu_exon/results/figures/coverage_summary/varLen1000/

smc_files="$( find ${input_res_dir} -name "summary_mean_cov.tsv" -type f | sort )"

# headers
echo -e "Min\t1st_Qu\tMedian\tMean\t3rd_Qu\tMax\tannot_list\tcond\trep\texp_des\tpath"

# results
for xxx in $smc_files; do
  samp_name="$( echo $xxx | awk -F '/' '{ print $(NF-2) }')"
  tail -n +2 $xxx | awk -v samp_name=$samp_name '{ OFS="\t"; print $0,samp_name }'
done


####
