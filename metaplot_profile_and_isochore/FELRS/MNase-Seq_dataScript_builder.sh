

script_fc () {
  local outdirsp=$1

  script_txt="
  public_script_dir=\"\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" && pwd )\"

  expSeq_bw_dir=\${data_dir}/bw_files/${link}
  source \${public_script_dir}/public_common_simpler.sh

  # maximum in heatmap
  ymax=100
  "

  echo "$script_txt"
}



# sed s%'lbdCtrl/'%''% | \
outdir=/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/data/bw_files
bwdir=~/analyses_SEBASTIEN/data/bigwig_files/public_MNase-Seq_MACS2_bw
scrdir=/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/coverage_bank

find ${bwdir}/ -name "*.bw" | \
sort | \
grep -vE 'lbdCtrl' | \
awk -F '/' '{ OFS="\t"; print $(NF-2),$NF }' | \
sed s%'_.*\t'%'\t'% | \
sed s%'\.bw$'%''% | \
while read line; do
  link="$( echo "$line" | awk '{ print "public_"$2"_GSE"$1 }' )"
  resfi="$( find ${bwdir}/$( echo $line | awk '{ print $1 }' )* -name "$( echo $line | awk '{ print $2 }' ).bw" )"
  if [[ "$( echo "$resfi" | wc -l )" != "1" ]]; then
    continue
  fi
  # echo "$link"
  # echo "$line"
  # echo "$resfi"

  outdirsp=${outdir}/$link
  mkdir -vp $outdirsp >&2
  ln -s $resfi ${outdirsp}/${link}.bw

  script_fc $outdirsp > ${scrdir}/${link}_script.sh
  # break
done
# | \
less


####
