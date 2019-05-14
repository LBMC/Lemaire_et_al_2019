#!/bin/bash

# Usage: bash main_total.sh [ --lists ] [ --plots ]

date

main_wdir=${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07

anal_dir="${main_wdir}/gene_metaVisu"
geneReg_dir=${anal_dir}
all_res_dir=${anal_dir}/res

for cline in ''; do # 'K562' 'HepG2' '293T' 'HeLa' 
  marq="_${cline}"
  if [ -z "$cline" ]; then
    marq=""
  fi

  exon_list_dir=${main_wdir}/gene_metaVisu/input_lists
  exon_list_pref=AT_rich${marq},GC_rich${marq}
  exon_lists=$( echo $exon_list_pref | sed s%'[0-9a-zA-Z_-]*'%"${exon_list_dir}/&\.txt"%g )

  annot_type_list=$( echo exon gene TSS_1kb geneEnd_1kb exon_first exon_last exon_intern intron | tr ' ' ',' ) #exon_ext  exon_intern_simple
  annot_type_list_total=exon_sync,${annot_type_list} #,exon_ssWinSmall
  
  # annot_type_list=intron_nm5 #intron_np0,intron_nm1,intron_nm2,intron_nm3,intron_nm4,intron_nm5,intron_np1,intron_np2,intron_np3,intron_np4 #'exon,exon_nm1,exon_nm2,exon_nm3,exon_nm4,exon_np1,exon_np2,exon_np3,exon_np4'
  # annot_type_list_total=intron_nm5_sync #exon_sync,exon_nm1_sync,exon_nm2_sync,exon_nm3_sync,exon_nm4_sync,exon_np1_sync,exon_np2_sync,exon_np3_sync,exon_np4_sync #intron_np0_sync,intron_nm1_sync,intron_nm2_sync,intron_nm3_sync,intron_nm4_sync,intron_nm5_sync,intron_np1_sync,intron_np2_sync,intron_np3_sync,intron_np4_sync #exon_sync,exon_intern,intron #exon,exon_nm1,exon_nm2,exon_nm3,exon_nm4,exon_np1,exon_np2,exon_np3,exon_np4

  col_pal='006600,3361FF'

  csvs_dir=${anal_dir}/csvs
  csvs="
GCp100_hg19_tracks.txt

# fig_13-01-19${marq}_tracks.txt
# fig_13-01-19${marq}_hg38_tracks.txt

# ChromatinLore_JBC_SF_ordered_list${marq}_hg19_missing_tracks.txt
# ChromatinLore_JBC_SF_ordered_list${marq}_hg38_missing_tracks.txt

# ChromatinLore_JBC_SF_ordered_list${marq}_hg19_tracks.txt
# ChromatinLore_JBC_SF_ordered_list${marq}_hg38_tracks.txt
# ChromatinLore_JBC_SF_ordered_list${marq}_hg18_tracks.txt

# Pol2_nucleosome_tracks_didier_26-11-18_hg19${marq}_tracks.txt
# Pol2_nucleosome_tracks_didier_26-11-18_hg38${marq}_tracks.txt

# punctual_tracks.txt
# 20-11-18${marq}_tracks.txt
# 20-11-18_hg38${marq}_tracks.txt
# 20-11-18_hg18${marq}_tracks.txt

# cistrome_Pol2_20-11-18${marq}_tracks.txt
# cistrome_20-11-18_tracks.txt
# Nuc_20-11-18${marq}_tracks.txt
# run_on${marq}_tracks.txt
# run_on_hg18${marq}_tracks.txt
# ChIP_Pol2_15-11-18${marq}_tracks.txt
# cline_POLR2${marq}_tracks.txt

# chromlore${marq}_tracks.txt
# cistrome_Pol2${marq}_tracks.txt
# cistrome${marq}_tracks.txt
# good_ChIP_H3_tracks.txt
# MNase_tracks.txt
# histoneMarks${marq}_tracks.txt
"
  csvs="$( echo $( echo "$csvs" | grep -vP "^(#|$)" | sed s%'^.*$'%"${csvs_dir}/&"% ) | tr ' ' ',' )"
  echo "$csvs"

  # selrs_parts_flag='--lists,--coverage'
  out_geneReg_dir=${geneReg_dir}/geneReg
  out_geneReg_tab_dir=${geneReg_dir}/geneReg_felrs
  out_geneReg_hg38_dir=${geneReg_dir}/geneReg_hg38
  out_geneReg_hg38_unlifted_dir=${geneReg_dir}/geneReg_hg38_unlifted
  out_geneReg_hg38_tab_dir=${geneReg_dir}/geneReg_hg38_felrs
  out_geneReg_hg18_dir=${geneReg_dir}/geneReg_hg18
  out_geneReg_hg18_unlifted_dir=${geneReg_dir}/geneReg_hg18_unlifted
  out_geneReg_hg18_tab_dir=${geneReg_dir}/geneReg_hg18_felrs

  mkdir -vp $out_geneReg_dir $out_geneReg_tab_dir $out_geneReg_hg38_dir $out_geneReg_hg38_tab_dir $out_geneReg_hg38_unlifted_dir $out_geneReg_hg18_dir $out_geneReg_hg18_tab_dir $out_geneReg_hg18_unlifted_dir >&2

  anal_name="geneVisu${marq}"

  if [[ "x${@}" == x*"--lists"* ]]; then
    ## prepare lists of genome regions for metaplots building
    for exon_list in $( echo $exon_lists | tr ',' ' ' ); do
      cmd="bash ../src/gene_metaVisu/filter_reg_list.sh $exon_list $out_geneReg_dir"
      echo "$cmd"
      eval "$cmd"

      for annot_type in $( echo ${annot_type_list} | tr ',' ' ' ); do
        for bed_file in $( ls ${out_geneReg_dir}/$( basename ${exon_list} .txt )_${annot_type}.bed ); do
          bed_hg38_file=${out_geneReg_hg38_dir}/$( basename $bed_file .bed ).hg38.bed
          out_unlift=${out_geneReg_hg38_unlifted_dir}/$( basename $bed_file .bed ).hg38_unlift.bed
          ../src/Liftover/liftOver <( sed s%'^[0-9XY]'%'chr&'% $bed_file ) ../src/Liftover/hg19ToHg38.over.chain.gz ${bed_hg38_file} ${out_unlift}
          sed -i '/^[^\t]*_\(alt\|random\)\t/d' ${bed_hg38_file}
          sed -i s%'^chr'%''% ${bed_hg38_file}
          sed -i s%'-$'%'-1'% ${bed_hg38_file}
          sed -i s%'+$'%'1'% ${bed_hg38_file}

          bed_hg18_file=${out_geneReg_hg18_dir}/$( basename $bed_file .bed ).hg18.bed
          out_unlift=${out_geneReg_hg18_unlifted_dir}/$( basename $bed_file .bed ).hg18_unlift.bed
          ../src/Liftover/liftOver <( sed s%'^[0-9XY]'%'chr&'% $bed_file ) ../src/Liftover/hg19ToHg18.over.chain.gz ${bed_hg18_file} ${out_unlift}
          sed -i '/^[^\t]*_\(alt\|random\)\t/d' ${bed_hg18_file}
          sed -i s%'^chr'%''% ${bed_hg18_file}
          sed -i s%'-$'%'-1'% ${bed_hg18_file}
          sed -i s%'+$'%'1'% ${bed_hg18_file}

          out_file=${out_geneReg_tab_dir}/$( basename $bed_file .bed ).tab
          echo -e "exon_id\tcoord\tstrand" > ${out_file}
          awk -F '\t' '{ OFS="\t"; print $4,$1":"$2+1"-"$3,$6 }' ${bed_file} >> ${out_file}

          out_hg38_file=${out_geneReg_hg38_tab_dir}/$( basename $bed_file .bed ).tab
          echo -e "exon_id\tcoord\tstrand" > ${out_hg38_file}
          awk -F '\t' '{ OFS="\t"; print $4,$1":"$2+1"-"$3,$6 }' ${bed_hg38_file} >> ${out_hg38_file}

          out_hg18_file=${out_geneReg_hg18_tab_dir}/$( basename $bed_file .bed ).tab
          echo -e "exon_id\tcoord\tstrand" > ${out_hg18_file}
          awk -F '\t' '{ OFS="\t"; print $4,$1":"$2+1"-"$3,$6 }' ${bed_hg18_file} >> ${out_hg18_file}
        done
      done
    done

    selrs_parts_flag='--lists-coverage'

    echo "bash ../src/gene_metaVisu/annot_lists_covSum.sh ${exon_list_pref} ${col_pal} ${csvs} ${annot_type_list_total} ${anal_name} ${selrs_parts_flag}"
    bash ../src/gene_metaVisu/annot_lists_covSum.sh ${exon_list_pref} ${col_pal} ${csvs} ${annot_type_list_total} ${anal_name} ${selrs_parts_flag} ${geneReg_dir} ${all_res_dir}
    bash ../src/gene_metaVisu/annot_lists_covSum.sh ${exon_list_pref} ${col_pal} ${csvs} ${annot_type_list_total} ${anal_name} ${selrs_parts_flag} ${geneReg_dir} ${all_res_dir} --hg38
    bash ../src/gene_metaVisu/annot_lists_covSum.sh ${exon_list_pref} ${col_pal} ${csvs} ${annot_type_list_total} ${anal_name} ${selrs_parts_flag} ${geneReg_dir} ${all_res_dir} --hg18
  fi
  # break


  if [[ "x${@}" == x*"--plots"* ]]; then
    ## launch coverage metaplots building
    selrs_parts_flag='--coverage'

    for csvs_file in $( echo $csvs | tr ',' ' ' ); do
      if [ -z "$( ls ${csvs_file} 2> /dev/null )" ]; then
        echo '!!! Does not exist: '"${csvs_file}" >&2
        continue
      fi

      hg38_flag=''
      if [[ "${csvs_file}" == *"cistrome"* || "${csvs_file}" == *"_hg38_"* ]]; then
        hg38_flag='--hg38'
      elif [[ "${csvs_file}" == *"_hg18_"* ]]; then
        hg38_flag='--hg18'
      fi
      cmd="bash ../src/gene_metaVisu/annot_lists_covSum.sh ${exon_list_pref} ${col_pal} ${csvs_file} ${annot_type_list_total} ${anal_name} ${selrs_parts_flag} ${geneReg_dir} ${all_res_dir} ${hg38_flag}"
      echo "$cmd"
      eval "$cmd"
    done
  fi

  # break
done


date
####
