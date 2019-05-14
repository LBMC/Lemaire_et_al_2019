
data_dir=../src/gene_metaVisu/data

sel_annot_onGeneId () {
  query=$1
  store_bed=$2

  join -t "$( printf '\t' )" -1 7 -2 1 <( paste ${store_bed} <( awk -F '\t' '{ OFS="\t"; print $4 }' $store_bed ) | sed s%'_[0-9]*$'%''% | sort -k7,7 ) <( less ${query} | awk '{ print $1 }' | sed s%'_[0-9]*'%''% | sort -k1,1 ) | \
  cut -d "$( printf '\t' )" -f 1 --complement | awk '!a[$0]++'
}

sel_annot_onExonId () {
  query_arg=$1
  store_bed_arg=$2

  join -t "$( printf '\t' )" -1 7 -2 1 <( paste ${store_bed_arg} <( awk -F '\t' '{ OFS="\t"; print $4 }' $store_bed ) | sort -k7,7 ) <( less ${query_arg} | awk '{ print $1 }' | sort -k1,1 ) | \
  cut -d "$( printf '\t' )" -f 1 --complement | awk '!a[$0]++'
}

nbh_compute () {
  query_arg=$1
  store_bed_arg=$2
  nbh_sign_arg=$3

  nsc="$( echo $nbh_sign_arg | awk -F '' '{ if ( $2 == "m" ) { $2="-" } else if ( $2 == "p" ) { $2="+" }; OFS=""; print $2,$3 }' )"

  head -n 1 ${query_arg}
  join -t "$( printf '\t' )" -1 1 -2 4 <( tail -n +2 ${query_arg} | awk -F '_' -v nsc=${nsc} '{ OFS="_"; print $1,$2+nsc}' | sort ) <( sort -k4,4 ${store_bed_arg} ) | cut -d "$( printf '\t' )" -f 1
}

sel_annotGene_onGeneId () {
  query=$1
  store_bed=$2
  join -t "$( printf '\t' )" -1 7 -2 1 <( paste ${store_bed} <( awk -F '\t' '{ OFS="\t"; print $4 }' $store_bed ) | tail -n +3 | sort -k7,7 ) <( less ${query} | awk '{ print $1 }' | awk -F '_' '{ print $1 }' | sort -k1,1 ) | \
  cut -d "$( printf '\t' )" -f 1 --complement | awk '!a[$0]++' | sed s%'^\ *'%''%
}


query=$1 #~/work_temp/rna_folding_measures/result/exon_lists/SFs/down_regulated_exons_by_SRSF9_in_K562.txt
out_dir=$2 #'.'
echo "$query" >&2

## bed for input exons
echo "exon" >&2
store_bed=${data_dir}/exon.bed
out_bed=${out_dir}/$( basename ${query} .txt)_$( basename ${store_bed} )
sel_annot_onExonId ${query} ${store_bed} > ${out_bed}

## bed for neighbour input exons
echo "neighbour exon" >&2
store_bed=${data_dir}/exon_intern.bed
for nbh_sign in nm4 nm3 nm2 nm1 np1 np2 np3 np4; do
  echo "  exon_${nbh_sign}" >&2
  inter_out=$( dirname $query )/nbh; mkdir -vp ${inter_out}
  inter_list=${inter_out}/$( basename ${query} .txt )_${nbh_sign}.txt
  out_bed=${out_dir}/$( basename ${query} .txt)_$( basename ${store_bed} _intern.bed )_${nbh_sign}.bed
  nbh_compute ${query} ${store_bed} ${nbh_sign} > ${inter_list}
  sel_annot_onExonId ${inter_list} ${store_bed} > ${out_bed}
done

## bed for neighbour input exons
echo "downstream intron" >&2
store_bed=${data_dir}/intron.bed
intron_out=$( dirname $query )/intron; mkdir -vp ${intron_out}
out_intron=${intron_out}/$( basename ${query} .txt )_intron_np0.txt
out_bed=${out_dir}/$( basename ${out_intron} .txt).bed
echo "intron_gene_pos" > ${out_intron}
sel_annot_onExonId ${query} ${store_bed} | awk -F '\t' '{ print $4 }' >> ${out_intron}
sel_annot_onExonId ${out_intron} ${store_bed} > ${out_bed}

echo "neighbour intron" >&2
for nbh_sign in nm5 nm4 nm3 nm2 nm1 np1 np2 np3 np4; do
  echo "  intron_${nbh_sign}" >&2
  inter_list=${intron_out}/$( basename ${query} .txt )_intron_${nbh_sign}.txt
  out_bed=${out_dir}/$( basename ${query} .txt)_intron_${nbh_sign}.bed
  nbh_compute ${query} ${store_bed} ${nbh_sign} > ${inter_list}
  sel_annot_onExonId ${inter_list} ${store_bed} > ${out_bed}
done

## bed for input genes
echo "gene" >&2
store_bed=${data_dir}/gene.bed
out_bed=${out_dir}/$( basename ${query} .txt)_$( basename ${store_bed} )
sel_annotGene_onGeneId ${query} ${store_bed} > ${out_bed}

## from genes containing query exons
store_name_list=( TSS_1kb exon_first exon_intern intron exon_last geneEnd_1kb )
for store_bed_name in ${store_name_list[*]}; do
  echo "$store_bed_name" >&2
  store_bed=${data_dir}/${store_bed_name}.bed
  out_bed=${out_dir}/$( basename ${query} .txt )_$( basename ${store_bed} )
  sel_annot_onGeneId ${query} ${store_bed} > ${out_bed}
done


####
