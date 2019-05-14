
results_dir=$1 #~ /scratch/slemaire/RNA-Seq_datasets/GSE74036_MCF7_epiTherapy_1/RNA-Seq/SINGLE
archive_dir=$2 #~ /Xnfs/site/auboeuf/data/SEBASTIEN/RNA-Seq_chromatinMarks # replace /scratch/slemaire/RNA-Seq_datasets
sign=$3 #~ pattern used to identify the files to store
# replace_dir=$3 #~ /scratch/slemaire/RNA-Seq_datasets

echo "dir for search: ${results_dir}"
echo "dir for storing: ${archive_dir}"
echo "signature to identify the file: ${sign}"

#~ find ${project_dir}/* -type f -o -type l > chemins_fich.txt
chemins_fich="$( find ${results_dir}/* -type f -o -type l )"
# echo "$chemins_fich" >&2

# patterns:
patterns="
/SS_sequences/.*${sign}.*
/CN_filtered/.*${sign}
/exons_lists/.*${sign}.*
/feature_tabs/${sign}.*_feat.tab
/splicing_site_regions/.*${sign}.*
/figures/coverage_summary/inEx[0-9]+_outOff[0-9]+/[^/]*/cov/[0-9]+clusters/[0-9]+mms/clustering_data/.*${sign}.*
"

for xxx in $patterns; do
  echo "   >>> $xxx" >&2
  for yyy in $( echo "$chemins_fich" | grep -E "$xxx" | sort ); do
    # echo "$yyy" # >&2
    new=$( echo $yyy | sed s,"${results_dir}","${archive_dir}", )
    mkdir -vp $( dirname $new )
    cp -a $yyy $new
  # echo ''
  done
done # | less
