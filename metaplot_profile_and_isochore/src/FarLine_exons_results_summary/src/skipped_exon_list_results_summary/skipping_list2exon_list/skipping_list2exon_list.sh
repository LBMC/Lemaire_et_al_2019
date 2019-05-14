#!/bin/bash

skipping_list2exon_list_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${skipping_list2exon_list_dir}/exon_coord_extraction.sh
source ${skipping_list2exon_list_dir}/exons_list2BED6_script.sh
source ${skipping_list2exon_list_dir}/exon2SSregion.sh

out_list_dir=$1 #~ ./results/exons_lists; mkdir -vp ${ex_list_dir}
feat_fasterdb=$2 #~ ./data/all_faster_db_web.tsv
listes_exons=$3 #~ /home/sebastien/work_temp/id_card_temp/results/all_data_tab/exon_up_siDNMT3b-siGL2.tab
bed_exons_genomiques_bis=$4 #~ ${data_dir}/exons_genomiques_bis.tsv
CONF_FILE=$5
CODING_EXON_TAB=$6
feat_tab_dir=$7
# listes_refs=$7
# names_refs=$8

listes_refs=''; names_refs=''
flag='--ex-tab-ref'
if [[ "x${@}" == x*"${flag}"* ]]; then
  listes_refs="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 1 )"
  names_refs="$( echo ${@} | awk -F "${flag} " '{ print $2 }' | cut -d ' ' -f 2 )"
fi

input_type_arg=''
if [[ "x${@}" == x*"--coord-only"* ]]; then
  input_type_arg="--coord-only"
fi

exon2SSregion_cmd () {
  local reg_list=$1
  local SS=$2
  local in_exon=$3
  local out_off=$4
  local reg_ss_list=$5
  cmd="Rscript ${skipping_list2exon_list_dir}/exon2SSregion.r ${reg_list} ${SS} ${in_exon} ${out_off} ${bed_exons_genomiques_bis} ${input_type_arg} > ${reg_ss_list}"
  echo "$cmd"
}

source $CONF_FILE

listes_refs_array=( $( echo ${listes_refs} | tr ',' ' ' ) )
names_refs_array=( $( echo ${names_refs} | tr ',' ' ' ) )
unset new_listes_refs_array; declare -A new_listes_refs_array
unset out_listes_refs_array; declare -A out_listes_refs_array
for ((iii=0; iii<${#listes_refs_array[@]}; iii++)); do
  name=${names_refs_array[$iii]}
  liste=${listes_refs_array[$iii]}
  new_listes_refs_array[${name}]=${liste}
  out_listes_refs_array[${name}]=${out_list_dir}/${name}_list.tsv
done

#~ exons_list=${out_list_dir}/$( basename $liste_exons .tab )_list.tsv
#~ exons_list=( $( echo $listes_exons | tr ',' '\n' | sed s/'.tab'$// | awk -F '/' -v old=${out_list_dir} '{ print old"/"$NF"_list.tsv" }' ) )
#~_unset exons_list_array; declare -A exons_list_array
listes_exons_array=( $( echo ${listes_exons} | tr ',' ' ' ) )
for xxx in ${listes_exons_array[*]}; do
  exons_list_array=( ${exons_list_array[*]} ${out_list_dir}/$( basename $xxx .tab )_list.tsv )
done

reg_list_array=( ${exons_list_array[*]} )
if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
  reg_list_array=( ${out_listes_refs_array[*]} ${reg_list_array[*]} )
fi


################################
###### For base fraction metagene
echo -e "\n>>> lists for base fraction metagenes" >&2
#### work on exons
echo -e "\n>>> lists of exons" >&2
## build the complete table of coordinates in one file
if [[ "x${@}" == x*"--coord-only"* ]]; then
  # for reference region lists
  for name in ${names_ref_array[*]}; do
    ln -fs ${new_listes_refs_array[$name]} ${out_listes_refs_array[$name]}
  done

  # for query regions
  for ((xxx=0;xxx<${#exons_list_array[@]};xxx++)); do
    ln -fs ${listes_exons_array[$xxx]} ${exons_list_array[$xxx]}
  done
else
  # for reference exon lists
  if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    for name in ${names_refs_array[*]}; do
      feature_tab=${feat_tab_dir}/${name}_feat.tab
      Rscript ${skipping_list2exon_list_dir}/ref_list.r $feat_fasterdb ${new_listes_refs_array[$name]} > ${feature_tab}
      exon_coord_extraction "FasterDB" ${feature_tab} > ${out_listes_refs_array[$name]}
    done
  fi

  # for query skipping exons/events
  for ((xxx=0;xxx<${#exons_list_array[@]};xxx++)); do
    exon_coord_extraction "FarLine" ${listes_exons_array[$xxx]} $bed_exons_genomiques_bis > ${exons_list_array[$xxx]}
  done
fi


## convert the list of exons ( already containing exon information as coordinates) in bed files
echo -e "\n>>> bed files of lists of exons" >&2
# out_list_bed_dir=${out_list_dir}/bed6; mkdir -vp ${out_list_bed_dir}
for reg_list in ${reg_list_array[*]}; do
  exons_list2BED6_script $reg_list #~  > ${out_list_bed_dir}/$( basename ${reg_list} .tsv ).bed
done


#### work on splicing sites
if [[ "x$@" == x*"--base-fraction"* ]]; then
  echo -e "\n>>> lists of splicing site regions, associated BED files, and extended regions" >&2
  in_exon_array=( $( echo ${in_exon_list} | tr ',' ' ' ) )
  out_off_array=( $( echo ${out_off_list} | tr ',' ' ' ) )
  for ((iii=0; iii<${#in_exon_array[@]}; iii++)); do
    in_exon=${in_exon_array[$iii]}
    out_off=${out_off_array[$iii]}

    out_ssRegion_dir=${out_list_dir}/../splicing_site_regions/inEx${in_exon}_outOff${out_off}; mkdir -vp ${out_ssRegion_dir}
    for reg_list in ${reg_list_array[*]}; do
      for SS in 3SS 5SS; do
        # build table of splicing site region coordinates
        reg_ss_list=${out_ssRegion_dir}/$( basename ${reg_list} _list.tsv )_${SS}.tsv
        # exon2SSregion ${reg_list} ${SS} ${in_exon} ${out_off} > ${reg_ss_list}
        # cmd="Rscript ${skipping_list2exon_list_dir}/exon2SSregion.r ${reg_list} ${SS} ${in_exon} ${out_off} ${bed_exons_genomiques_bis} ${input_type_arg} > ${reg_ss_list}"
        cmd="$( exon2SSregion_cmd ${reg_list} ${SS} ${in_exon} ${out_off} ${reg_ss_list} )"
        # echo "$cmd"
        eval "$cmd"

        # make the corresponding bed files
        # out_ssRegion_bed_dir=${out_ssRegion_dir}/bed6; mkdir -vp ${out_ssRegion_bed_dir}
        exons_list2BED6_script ${reg_ss_list} #~ > ${out_ssRegion_bed_dir}/$( basename ${reg_ss_list} .tsv ).bed

        # extend the splicing site region for computing in a window centered on each base
        for window_size in $( echo ${window_sizes} | tr ',' ' ' ); do
          out_ssRegion_ext_dir=${out_ssRegion_dir}/win${window_size}b; mkdir -vp ${out_ssRegion_ext_dir}
          out_ssRegion_ext_file=${out_ssRegion_ext_dir}/$( basename ${reg_ss_list} .tsv )_ext.tsv
          cmd="Rscript ${skipping_list2exon_list_dir}/reg_ss_coord.r ${reg_ss_list} ${window_size} > ${out_ssRegion_ext_file}"
          # echo "$cmd"
          eval "$cmd"
          exons_list2BED6_script ${out_ssRegion_ext_file}
        done
      done
    done
  done
fi

################################
###### For coverage metagene
if [[ "x$@" == x*"--coverage"* ]]; then
  echo -e "\n>>> lists for coverage and RRBS metagenes" >&2
  out_ssRegion_dir=${out_list_dir}/../splicing_site_regions/inEx${cov_in_exon}_outOff${cov_out_off}; mkdir -vp ${out_ssRegion_dir}
  for reg_list in ${reg_list_array[*]}; do
    for SS in 3SS 5SS; do
      reg_ss_list_cov=${out_ssRegion_dir}/$( basename ${reg_list} _list.tsv )_${SS}.tsv
      # exon2SSregion ${reg_list} ${SS} ${cov_in_exon} ${cov_out_off} > ${reg_ss_list_cov}
      Rscript ${skipping_list2exon_list_dir}/exon2SSregion.r ${reg_list} ${SS} ${cov_in_exon} ${cov_out_off} ${bed_exons_genomiques_bis} ${input_type_arg} > ${reg_ss_list_cov}
      cmd="$( exon2SSregion_cmd ${reg_list} ${SS} ${cov_in_exon} ${cov_out_off} ${reg_ss_list_cov} )"
      # echo "$cmd"
      eval "$cmd"
      exons_list2BED6_script ${reg_ss_list_cov}

      for window_size in $( echo ${RRBS_windows} | tr ',' ' ' ); do
        out_ssRegion_ext_dir=${out_ssRegion_dir}/win${window_size}b; mkdir -vp ${out_ssRegion_ext_dir}
        out_ssRegion_ext_file=${out_ssRegion_ext_dir}/$( basename ${reg_ss_list_cov} .tsv )_ext.tsv
        Rscript ${skipping_list2exon_list_dir}/reg_ss_coord.r ${reg_ss_list_cov} ${window_size} > ${out_ssRegion_ext_file}
        exons_list2BED6_script ${out_ssRegion_ext_file}
      done
    done

  done

  out_ssRegion_dir=${out_list_dir}/../splicing_site_regions/aroundCenter${cov_center_off}; mkdir -vp ${out_ssRegion_dir}
  for reg_list in ${reg_list_array[*]}; do
    reg_ss_list_cov=${out_ssRegion_dir}/$( basename ${reg_list} _list.tsv )_center.tsv
    cmd="Rscript ${skipping_list2exon_list_dir}/exonCenterRegion.r ${reg_list} ${cov_center_off} ${bed_exons_genomiques_bis} ${input_type_arg} > ${reg_ss_list_cov}"
    # echo "$cmd"
    eval "$cmd"

    exons_list2BED6_script ${reg_ss_list_cov}

    for window_size in $( echo ${RRBS_windows} | tr ',' ' ' ); do
      out_ssRegion_ext_dir=${out_ssRegion_dir}/win${window_size}b; mkdir -vp ${out_ssRegion_ext_dir}
      out_ssRegion_ext_file=${out_ssRegion_ext_dir}/$( basename ${reg_ss_list_cov} .tsv )_ext.tsv
      Rscript ${skipping_list2exon_list_dir}/reg_ss_coord.r ${reg_ss_list_cov} ${window_size} > ${out_ssRegion_ext_file}
      exons_list2BED6_script ${out_ssRegion_ext_file}
    done

  done
fi
