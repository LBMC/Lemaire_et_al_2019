
report_img_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

wdir=${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/res_jbc_anyCell
list_dir=${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/gene_metaVisu/report_cov_lists

signature='anyCell' #'byCell'

out_pertrack=./img_perTrack
mkdir -vp ${out_pertrack} >&2

out_final=./reports
mkdir -vp ${out_final} >&2


signal_line_sel () {
  fich=$1
  grep -vP '^($|#)' ${fich}
}

unset visu_type; declare -A visu_type
visu_type['TSS_1kb']=aroundCenter1000
visu_type['exon']=aroundCenter100
visu_type['exon_sync']=varLen1000
visu_type['exon_first']=varLen1000
visu_type['exon_intern']=varLen1000
# visu_type['exon_intern_simple']=varLen1000
visu_type['intron']=varLen1000
visu_type['exon_last']=varLen1000
visu_type['geneEnd_1kb']=aroundCenter1000
visu_type['gene']=varLen2000


# annot_type_list=( TSS_1kb exon_first exon_intern intron exon_last geneEnd_1kb )
annot_type_list=( TSS_1kb exon_first exon exon_sync exon_intern intron exon_last geneEnd_1kb gene )
# annot_type_list=( TSS_1kb exon_first exon_intern intron exon_last geneEnd_1kb gene )

signal_group_list='
# ChromatinLore_JBC_SF_hgxx.filtered.ordered_20181210.txt
ChromatinLore_JBC_SF_hgxx.20181210.txt
'
signal_group_list="$( echo "${signal_group_list}" | grep -vP '^(#|$)' | tr '\n' ' ' )"

measure_type='mean'

colors=( '66FF66' 'FF6666' )
col_nbr=0
# target_mem=''
sf_list="$( ls geneReg_jbc/*_down_exon.bed | awk -F '/' '{ print $NF }' | awk -F '_down' '{ print $1 }' | sed s%'_\(HEK293\|HepG2\|K562\|MCF7\|HeLa\)'%''% | sort | uniq | tr '\n' ' ' | less )"
work=0
for mysf in SNRNPC; do #${sf_list}; do
  echo "mysf: $mysf"

  for target in H2A H3 H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K4me3B H3K79me2 H3K79me3 H3K9ac H3K9me2 H3K9me3 H4K20me1 MNAse POLR2A POLR2AphosphoS2 POLR2AphosphoS5; do

    report_signature="${mysf}_${target}_${signature}" #_anyCell"
    if [ -n "$report_signature" ]; then
      report_signature="${report_signature}_"
    fi

    #if [[ "${work}" == "0" ]]; then
    #  if [[ "$mysf" == "HNRNPU" ]]; then
    #    work=1
    #  fi
    #  continue
    #fi

    for signal_group in ${signal_group_list}; do #histoneAc_clineChromLore histoneMeth_complete_refined
      echo -e "\n>>> $signal_group" >&2
      nbr_sig=$( signal_line_sel ${list_dir}/${signal_group} | wc -l )

      res_dir=${wdir}
      if [[ "x$signal_group" == x"cistrome"* ]]; then
        res_dir=${wdir}_hg38
      fi

      # cmd_loc="echo \${signal_list_${signal_group}[*]} | tr ',' ' '"
      # nbr_sig="$( eval $cmd_loc | tr ' ' '\n' | wc -l )"
      unset out_img_array; declare -a out_img_array
      nnn=0
      while read signal; do
        #~ echo "!!! input SIGNAL: $signal"
        ((nnn++))
        #echo -ne "$( python3 -c "print( int( ( ${nnn}/${nbr_sig} )*100 ) )" )/100\r"
        unset img_array; declare -a img_array

        for annot_type in ${annot_type_list[*]}; do #TSS_1kb exon_first exon_intern intron exon_last geneEnd_1kb
        vt_key=${annot_type}
          if [[ "$annot_type" == 'exon_sync' ]]; then
            annot_type='exon'
          fi
          # echo "--- $annot_type"
          # echo "${visu_type[${annot_type}]}"
          cmd="find ${wdir}/${mysf}_jbc*/geneVisu_${annot_type}/results/figures/coverage_summary/${visu_type[${vt_key}]}/${signal}/cov/metagene_${measure_type}/ -name \"metagene_${measure_type}*.png\" | grep _${target}_ "
          # echo ">>> $cmd"
          img="$( eval "$cmd" 2> /dev/null )"
          #~ echo -e "images:\n$img"
          if [ -z "$img" ]; then
            img="${report_img_dir}/blank_img_480x480p.png"
          fi

          if [[ "$( echo $img | wc -l )" -gt 1 ]]; then
            echo "2 or more images detected for `$signal`" >&2
          elif [[ -z $( ls -d ${wdir}/${mysf}_jbc*/geneVisu_${annot_type}/results/figures/coverage_summary/${visu_type[${vt_key}]}/${signal}/cov/metagene_${measure_type} 2> /dev/null ) && -n $( ls ${wdir}/${mysf}_jbc*/geneVisu_${annot_type}/results/figures/coverage_summary/${visu_type[${vt_key}]}/${signal}/cov/qsub_script.sh 2> /dev/null ) ]]; then
            echo -e "\t »»» qsub $( ls ${wdir}/${mysf}_jbc*/geneVisu_${annot_type}/results/figures/coverage_summary/${visu_type[${vt_key}]}/${signal}/cov/qsub_script.sh )"
          fi

          img_array=( ${img_array[*]} $img )
        done

        if [[ -z "$( echo "${img_array[*]}" | sed s%'[^\ ]*blank_img[^\ ]*'%''%g | sed s%'[\ ]*'%''%g )" ]]; then
          continue
        fi

        # target="$( echo "$signal" | awk -F '_' '{ print $3 }' )"
        # if [[ "$signal_group" == "histoneMeth" ]]; then
        #   target="$( echo "$signal" | sed s%'^\(GSE[0-9]*\|ENC[A-Z0-9]*\)_\(ENC[A-Z0-9]*_\)\{0,1\}\([12]_\|igv_\|noCTRL_\)\{0,1\}'%''% | awk -F '_' '{ print $2 }' )"
        # fi
        # echo $target >&2

        # if [[ "$target" != "$target_mem" ]]; then
        #   target_mem="$target"
        #   col_nbr=$((1-$col_nbr))
        # fi

        # echo "!!! SIGNAL: $signal"
        out_img=${out_pertrack}/${report_signature}${measure_type}_${signal}${marq}.png
        out_img_array=( ${out_img_array[*]} ${out_img} )
        cmd="montage -tile ${#annot_type_list[*]}x1 -geometry +0+0 ${img_array[*]} ${out_img}
        montage -label ${signal} -pointsize 30 -background '#${colors[$col_nbr]}' -geometry +0+0 -tile 1x1 ${out_img} ${out_img}"
        #~ echo "$cmd" | tr ' ' '\n'
        eval "$cmd"

      done < <( signal_line_sel ${list_dir}/${signal_group} )

      out_report=${out_final}/report_signal_${report_signature}${measure_type}_${signal_group}${marq}.pdf
      if [ "${#out_img_array[*]}" -lt 60 ]; then
        cmd="montage -tile 1x6 -geometry +0+0 ${out_img_array[*]} ${out_report}"
        echo "$cmd" | tr ' ' '\n'
        eval "$cmd"
      else
        nnn=0
        out_img_list=''
        out_img_idx=0
        for idx in $( seq 0 $((${#out_img_array[*]}-1)) ); do
          out_img_list="${out_img_list} ${out_img_array[$idx]}"
          ((nnn++))
          
          if [[ "$nnn" == 6 || "$((idx+1))" == "${#out_img_array[*]}" ]]; then
            #~ out_img_part=${out_pertrack}/$( basename ${out_report} )_${out_img_idx}
            out_img_part=$( dirname ${out_report} )/$( basename ${out_report} .pdf )_${out_img_idx}.pdf
            #~ out_img_part_array=( ${out_img_part_array[*]} ${out_img_part} )
            cmd="montage -tile 1x6 -geometry +0+0 ${out_img_list} ${out_img_part}"
            echo "$cmd" | tr ' ' '\n'
            eval "$cmd"
            ((out_img_idx++))

            nnn=0
            out_img_list=''
          fi
        done
      fi
    done
  done
done




####
