#!/bin/bash

## compute distribution of the exons in the isochores/TADs/…
Rscript ../src/exon_in_isochore/nbr_exon_in_TADs.r


#############################################
#############################################
## compute the %GC in the exons, isochrores, TADs, …
fasta_dir=../src/FarLine_exons_results_summary/data/hg19/EnsEMBL_all
tad_bed_dir=/home/sebastien/analyses_SEBASTIEN/data/misc/genome_organization/TADs
isochore_bed_dir=/home/sebastien/analyses_SEBASTIEN/data/misc/genome_organization/ISOCHORES
lad_bed_dir=/home/sebastien/analyses_SEBASTIEN/data/misc/genome_organization/LADs
heterochromatin_bed_dir=/home/sebastien/analyses_SEBASTIEN/data/misc/genome_organization/HETEROCHROMATIN
annot_bed_dir=/home/sebastien/analyses_SEBASTIEN/data
# bed_dir=${bed_dir}/short_bed
annot_file_list=$( echo $( find ${tad_bed_dir}/* -maxdepth 0 -name "*.bed" -type f ) | tr ' ' ',' )
annot_file_list="${annot_file_list},$( echo $( find ${isochore_bed_dir}/* -maxdepth 0 -name "*.bed" -type f ) | tr ' ' ',' )"
annot_file_list="${annot_file_list},$( echo $( find ${heterochromatin_bed_dir}/* -maxdepth 0 -name "*tsv.bed" -type f ) | tr ' ' ',' )"
annot_file_list="${annot_file_list},$( echo $( find ${lad_bed_dir}/* -maxdepth 0 -name "*.bed" -type f | grep -P '(Activated_Jurkat|Resting_Jurkat|lung_fibroblasts_lifted)' ) | tr ' ' ',' )"
annot_file_list="${annot_file_list},${annot_bed_dir}/exon.bed,${annot_bed_dir}/gene.bed"
# annot_file_list="${isochore_bed_dir}/ISOFINDER_isochore_hg19_reComp_methP3.bed"
# annot_file_list=$( echo $( echo $( echo $annot_file_list | tr ',' '\n' | tail -2 ) | tr '\n' ' ' ) | tr ' ' ',' )
out_dir=results/GC_content
mkdir -vp $out_dir
out_file_list=$( echo $annot_file_list | sed s%"${tad_bed_dir}"%"${out_dir}"%g | sed s%"${isochore_bed_dir}"%"${out_dir}"%g | sed s%"${heterochromatin_bed_dir}"%"${out_dir}"%g | sed s%"${lad_bed_dir}/[a-zA-Z0-9]*"%"${out_dir}"%g | sed s%"${annot_bed_dir}"%"${out_dir}"%g | sed s%'\.bed'%'\.seq'%g )

cust_fc_file=../src/exon_in_isochore/fract_gc.py

# python3 ../src/exon_in_isochore/list_coord_seq_retriever.py $fasta_dir $annot_file_list $out_file_list \
#   --bed6 --remove-chr --cust-fc $cust_fc_file

python3 ../src/exon_in_isochore/list_coord_seq_retriever.py $fasta_dir $annot_file_list $out_file_list \
  --bed6 --remove-chr --without-header --cust-fc $cust_fc_file --nbr-cores 4


## extract AT-rich and GC-rich exons/genes from the whole pool
gcat_seq=${out_dir}/GC-AT-rich_exon.seq
exon_seq=${out_dir}/exon.seq
head -1 $exon_seq > $gcat_seq
gcat_txt_files=( ../gene_metaVisu/input_lists/GC_rich.txt ../gene_metaVisu/input_lists/AT_rich.txt ./CCE.txt )
for base_rich in ${gcat_txt_files[*]}; do
  out_seq_file=${out_dir}/$( basename $base_rich .txt )_exon.seq
  head -1 $exon_seq > ${out_seq_file}
  join -t "$( printf '\t' )" -1 1 -2 4 <( sort $base_rich ) <( sort -k4,4 $exon_seq ) | awk '{ OFS="\t"; print $2,$3,$4,$1 }' >> ${out_dir}/$( basename $base_rich .txt )_exon.seq

  if [[ "$( basename $base_rich )" == 'GC_rich.txt'* || "$( basename $base_rich )" == 'AT_rich.txt'* ]]; then
    tail -n +2 ${out_seq_file} >> $gcat_seq
  fi
done


exonid_to_geneid_fc () {
  gcat_seq=$1
  less $gcat_seq | tail -n +2 | awk '{ print $NF }' | awk -F '_' '{ print $1 }' | sort -n | uniq
}

gcat_gene=${out_dir}/GC-AT-rich_exon_gene.txt
echo "exon_gene_pos" > $gcat_gene
exonid_to_geneid_fc $gcat_seq >> $gcat_gene


for base_rich in ${gcat_txt_files[*]}; do
  out_gene_id=${out_dir}/$( basename $base_rich .txt )_exon_gene.txt
  exon_seq_file=${out_dir}/$( basename $base_rich .txt )_exon.seq
  echo "exon_gene_pos" > $out_gene_id
  exonid_to_geneid_fc $exon_seq_file >> $out_gene_id
done

gene_seq_select_fc () {
  gcat_gene=$1
  gene_seq=$2
  join -t "$( printf '\t' )" -1 1 -2 4 <( sort $gcat_gene ) <( sort -k4,4 $gene_seq ) | awk '{ OFS="\t"; print $2,$3,$4,$1 }'
}

gcat_gene_seq=${out_dir}/GC-AT-rich_exon_gene.seq
gene_seq=${out_dir}/gene.seq
head -1 $gene_seq > $gcat_gene_seq
gene_seq_select_fc $gcat_gene $gene_seq >> $gcat_gene_seq
# join -t "$( printf '\t' )" -1 1 -2 4 <( sort $gcat_gene ) <( sort -k4,4 $gene_seq ) | awk '{ OFS="\t"; print $2,$3,$4,$1 }' >> $gcat_gene_seq

gcat_gene_seq=${out_dir}/GC_rich_exon_gene.seq
gene_seq=${out_dir}/gene.seq
head -1 $gene_seq > $gcat_gene_seq
gene_seq_select_fc ${out_dir}/GC_rich_exon_gene.txt $gene_seq >> $gcat_gene_seq

gcat_gene_seq=${out_dir}/AT_rich_exon_gene.seq
gene_seq=${out_dir}/gene.seq
head -1 $gene_seq > $gcat_gene_seq
gene_seq_select_fc ${out_dir}/AT_rich_exon_gene.txt $gene_seq >> $gcat_gene_seq

gcat_gene_seq=${out_dir}/CCE_exon_gene.seq
gene_seq=${out_dir}/gene.seq
head -1 $gene_seq > $gcat_gene_seq
gene_seq_select_fc ${out_dir}/CCE_exon_gene.txt $gene_seq >> $gcat_gene_seq


#############################################
#############################################
## convert in BED6 format
for gc_file in $( find ${out_dir}/* -maxdepth 0 -name "*.seq" ); do
  # gc_file=${out_dir}/$( basename $xxx .bed ).seq
  out_bed=${out_dir}/$( basename $gc_file .seq ).bed
  less $gc_file | awk '{ OFS="\t"; print $1,$4,$3,$2 }' | tr ':' '\t' | tr '-' '\t' | awk '{ OFS="\t"; if ( $6 == "1" ) { $6 = "+" }; if ( $6 == "-1" ) { $6 = "-" }; print $0 }'> $out_bed
done


#############################################
#############################################
## plots / corrélations / counts
#### function
GCcontent_tabs_to_2DplotDensity_fc () {
  tab1=$1
  tab2=$2
  plot_dir=$3
  bdt_int_dir=$4
  count_dir=$5

  base1=$( basename $tab1 .seq )
  base2=$( basename $tab2 .seq )

  out_sign=${base1}_vs_${base2}
  axes_arg="--axes-labs ${base1},${base2}"

  ## density plot and statistics
  Rscript ../src/exon_in_isochore/GCcontent_tabs_to_2DplotDensity.r $tab1 $tab2 $plot_dir $bdt_int_dir $out_sign --add-chr ${axes_arg}

  ## overlap counts
  #~ python3 ../src/exon_in_isochore/bdt_int_to_ovlCounts.py ${bdt_int_dir}/${out_sign}.tab > ${count_dir}/${out_sign}.txt
}

####

plot_dir=results/GC_content/2D_plots
bdt_int_dir=results/GC_content/intersections
count_dir=results/GC_content/counts
mkdir -vp $plot_dir $bdt_int_dir ${count_dir}


# exon/gene Vs isochores/TADs/…
tab1_list="
${out_dir}/GC-AT-rich_exon.seq"
# ${out_dir}/exon.seq
# ${out_dir}/gene.seq
# ${out_dir}/GC-AT-rich_exon_gene.seq
# ${out_dir}/GC_rich_exon_gene.seq
# ${out_dir}/GC_rich_exon.seq
# ${out_dir}/AT_rich_exon_gene.seq
# ${out_dir}/AT_rich_exon.seq
# ${out_dir}/CCE_exon_gene.seq
# ${out_dir}/CCE_exon.seq

tab2_list="$( find ${out_dir}/* -maxdepth 0 -name "*.seq" -type f | grep -vP '(exon|gene)' )"

for tab1 in $tab1_list; do
  echo ">>> $tab1"
  for tab2 in $tab2_list; do
    echo "--- $tab2"
    # tab1=results/GC_content/exon.seq
    # tab2=results/GC_content/GSE35156_RenLab_TADs_IMR90.lifted.hg19.seq

    cmd="GCcontent_tabs_to_2DplotDensity_fc $tab1 $tab2 $plot_dir $bdt_int_dir $count_dir"
    echo "$cmd"
    eval "$cmd"
    # break
  done
  # break
done

# exon Vs gene
for tab1 in ${out_dir}/exon.seq ${out_dir}/GC-AT-rich_exon.seq ${out_dir}/AT_rich_exon.seq ${out_dir}/GC_rich_exon.seq ${out_dir}/CCE_exon.seq; do
  # tab1=results/GC_content/exon.seq
  tab2=results/GC_content/gene.seq

  cmd="GCcontent_tabs_to_2DplotDensity_fc $tab1 $tab2 $plot_dir $bdt_int_dir $count_dir"
  echo "$cmd"
  eval "$cmd"
done



##################################
##################################
## TABLES for Didier / COUNTS of exons/genes in isochores/TADs/…
bash ../src/exon_in_isochore/counts_exonOrGene_in_TAD.sh

for xxx in $( find ${out_dir}/*.didier ); do
  less $xxx | awk '{ if ( $6 != "0" ) print }' > ${xxx}.short
done



##################################
##################################
## CORRELATION %GC of TAD with fraction of GC exon ( among AT/GC exons )
corr_out_dir=results/GC_content/corr_isoGCp100_exonGC
mkdir -vp $corr_out_dir

for min_nbr_exon_atgc in 0 5; do
  for iso_exon_counts in $( find results/GC_content/* -maxdepth 0 -name "*.didier" ); do
    Rscript ../src/exon_in_isochore/corr_isoGCp100_nbrExonGC.r ${iso_exon_counts} ${corr_out_dir} ${min_nbr_exon_atgc}
    # break
  done
  # break
done



##################################
##################################
## REPART exon AT/GC per isochore category
unset iso_array; declare -gA iso_array
iso_array['Costantini']=results/GC_content/ISO_Costantini.lifted.hg19.seq.didier
iso_array["Isosegmenter"]=results/GC_content/ISO_isosegmenter.hg19.seq.didier
iso_array["Isofinder"]=results/GC_content/ISOFINDER_isochores_hg19.seq.didier

out_repart_dir=results/GC_content/repart_plots
mkdir -vp $out_repart_dir

for min_nbr_exon_atgc in 0 5; do
  for name in ${!iso_array[*]}; do
    Rscript ../src/exon_in_isochore/repart_exon_isoCat.r ${iso_array[$name]} $out_repart_dir $name $min_nbr_exon_atgc
  done
done




####
