#!/bin/bash

base_fraction_metagene_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${base_fraction_metagene_dir}/splicingSiteSeq/DNA_seq_retriever_uniq.sh
source ${base_fraction_metagene_dir}/metagene_plots/metagene_plots_script.sh
source ${base_fraction_metagene_dir}/seqFile2Fasta.sh
source ${base_fraction_metagene_dir}/weblogo/weblogo_plots_script.sh
source ${base_fraction_metagene_dir}/rna_folding/rna_folding.sh
source ${base_fraction_metagene_dir}/methylation_metagene/methylation_extractor.sh
source ${base_fraction_metagene_dir}/methylation_metagene/methylation_metagene_script.sh


base_fraction_metagene () {
  local prefixes=$1 # siDNMT3a-siGL2_up,siDNMT3a-siGL2_down
  local FASTA_DIR=$2 # data/hg19/EnsEMBL_all
  local SS=$3 # 3SS
  local metagene_param=$4 # 20,50,100
  local prefixes_ref=$5 # psi_100_95_FDB,psi_60_40_FDB,psi_5_0_FDB
  local ssRegion_ext_dir=$6 # results_psi_FDB/splicing_site_regions/inEx50_outOff100/win1b
  local min_reads=$7 # minimal number of aligned reads on a CpG to considered it in mCpG_metagene

  local NBRCORES=1
  if [[ "x$@" == x*"--nbr-cores"* ]]; then
    NBRCORES=$( echo "$@" | awk -F '--nbr-cores' '{ print $NF }' | cut -d ' ' -f 2 )
  fi

  # check for a defined figure prefix
  local fig_prefix='test'
  if [[ "x${@}" == x*"--fig-prefix"* ]]; then
    local fig_prefix="$( echo ${@} | awk -F '--fig-prefix ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi

  # check for specific colors in metagene plots
  if [[ "x${@}" == x*"--color-pallette"* ]]; then
    local color_pallette="--color-pallette $( echo ${@} | awk -F '--color-pallette ' '{ print $2 }' | cut -d ' ' -f 1 )"
  fi


  ## define output_directories
  local fig_dir=$( echo ${ssRegion_ext_dir} | sed s,splicing_site_regions,figures/base_fraction, ); mkdir -vp ${fig_dir}
  local meth_meta_dir=$( echo ${fig_dir} | sed s,figures/base_fraction,figures/methCpG_meta, ); mkdir -vp ${meth_meta_dir}
  local out_seq_dir=$( echo ${ssRegion_ext_dir} | sed s,splicing_site_regions,SS_sequences, ); mkdir -vp ${out_seq_dir}
  local out_fa_dir=${out_seq_dir}/fasta; mkdir -vp ${out_fa_dir}
  local out_meth_tab_dir=${out_seq_dir}/CpG_meth; mkdir -vp ${out_meth_tab_dir}
  local out_rna_folding_dir=$( echo ${ssRegion_ext_dir} | sed s,splicing_site_regions,rna_folding, ); mkdir -vp ${rna_folding_dir}
  local rna_folding_fig_dir=$( echo ${fig_dir} | sed s,figures/base_fraction,figures/rna_folding, ); mkdir -vp ${rna_folding_fig_dir}

  ## recover the coordinates of extended splicing regions files
  local ssRegion_ext_list=${ssRegion_ext_dir}/$( echo ${prefixes} | sed s?','?"_${SS}_ext.tsv,${ssRegion_ext_dir}/"?g )_${SS}_ext.tsv
  local ssRegion_ext_ref_list=${ssRegion_ext_dir}/$( echo ${prefixes_ref} | sed s?','?"_${SS}_ext.tsv,${ssRegion_ext_dir}/"?g )_${SS}_ext.tsv

  # add the ref exon lists to the treatment
  if [[ "x${@}" == x*"--with-refs-exons"* ]]; then
    ssRegion_ext_list="${ssRegion_ext_ref_list},${ssRegion_ext_list}"
  fi


  # retrieve the sequences of the splicing regions for the query exons
  echo -e "\n>>> retrieve the sequences of splicing site regions: query list" >&2
  seq_file_array=( '' )
  for ssRegion_ext in $( echo ${ssRegion_ext_list} | tr ',' ' ' ); do
    seq_file_array=( ${seq_file_array[*]} ${out_seq_dir}/$( basename ${ssRegion_ext} .tsv ).seq )
  done
  seq_files=$( echo ${seq_file_array[*]} | tr ' ' ',' )
  cmd="
  DNA_seq_retriever_uniq $ssRegion_ext_list $FASTA_DIR $seq_files
  "
  echo "$cmd"
  # eval "$cmd"


  prefixes_list=${prefixes_ref},${prefixes}
  seq_files_all=()
  fasta_files_all=()
  for xxx in $( echo ${prefixes_list} | tr ',' ' ' ); do
    seq_file=${out_seq_dir}/${xxx}_${SS}_ext.seq
    seq_files_all=( ${seq_files_all[*]} ${seq_file} )

    fasta_file=${out_fa_dir}/$( basename ${seq_file} .seq ).fa
    fasta_files_all=( ${fasta_files_all[*]} ${fasta_file} )
  done

  # create the corresponding fasta file
  cmd="
  seqFile2Fasta $( echo ${seq_files_all[*]} | tr ' ' ',' ) $( echo ${fasta_files_all[*]} | tr ' ' ',' )
  "
  echo "$cmd"
  # eval "$cmd"


  #### the metagene of base fractions, of CpG methylation, the RNAfolding free energies, and the weblogos
  # set the lists of files

  ## build the metagene of nucleotide composition
  echo -e "\n>>> build the metagene of base fractions" >&2
  cmd="
  metagene_plots_script $( echo ${seq_files_all[*]} | tr ' ' ',' ) ${fig_dir} ${SS} ${metagene_param} ${prefixes_list} \
    --nbr-cores $NBRCORES --fig-prefix ${fig_prefix} ${color_pallette}
  "
  echo "$cmd"
  # eval "$cmd"

  ## build the weblogos of nucleotide composition
  echo -e "\n>>> build the weblogos" >&2
  cmd="
  weblogo_plots_script $( echo ${fasta_files_all[*]} | tr ' ' ',' ) ${fig_dir} ${SS} ${metagene_param} ${prefixes_list} \
    --nbr-cores $NBRCORES
  "
  echo "$cmd"
  # eval "$cmd"

  ## build the repartition plots of RNA folding free energies
  echo -e "\n>>> build repartition plots of RNA folding free energies" >&2
  cmd="
  rna_folding $( echo "${fasta_files_all[*]}" | tr ' ' ',' ) ${prefixes_list} ${out_rna_folding_dir} ${rna_folding_fig_dir} \
    --nbr-cores $NBRCORES --fig-prefix ${fig_prefix}
  "
  echo "$cmd"
  # eval "$cmd"


  #### CpG methylation
  echo -e "\n>>> retrieve the CpG methylation rates of splicing site regions: query list" >&2
  for expSeq_source in $( find ${base_fraction_metagene_dir}/methylation_metagene/coverage_bank/ -name "*_script.sh" | sort ); do # | grep WGBS ); do
    echo $( basename $expSeq_source ) >&2
    source $expSeq_source

    local meth_out_dir=${out_meth_tab_dir}/$( basename ${expSeq_source} _script.sh )
    mkdir -vp ${meth_out_dir} >&2

    cond_array=( $( echo $expSeq_cond | tr ',' ' ' ) )
    rep_array=( $( echo $expSeq_rep | tr ',' ' ' ) )
    bam_array=( $( echo $expSeq_bam_files | tr ',' ' ' ) )
    for (( spl_id=0; spl_id<${#bam_array[*]}; spl_id++ )); do
      # extract the CpG methylations
      cond=${cond_array[$spl_id]}
      rep=${rep_array[$spl_id]}
      bam_file=${bam_array[$spl_id]}

      meth_file_array=()
      for seq_file in ${seq_file_array[*]}; do
        meth_file_array=( ${meth_file_array[*]} ${meth_out_dir}/$( basename ${seq_file} _ext.seq )_${cond}n${rep}.tsv )
      done

      cmd="
      methylation_extractor ${bam_file} $( echo ${seq_file_array[*]} | tr ' ' ',' ) $( echo ${meth_file_array[*]} | tr ' ' ',' )
      "
      echo "$cmd"
      # eval "$cmd"
    done


    # build the methylation metagenes
    echo -e "\n>>> build the metagenes of CpG methylation" >&2
    local out_meth_meta_dir=${meth_meta_dir}/$( basename ${expSeq_source} _script.sh )
    mkdir -vp ${out_meth_meta_dir} >&2

    all_meth_file_array=()
    for (( spl_id=0; spl_id<${#bam_array[*]}; spl_id++ )); do
      # extract the CpG methylations
      cond=${cond_array[$spl_id]}
      rep=${rep_array[$spl_id]}
      bam_file=${bam_array[$spl_id]}

      meth_file_array=()
      for seq_file in ${seq_files_all[*]}; do
        local meth_file=${meth_out_dir}/$( basename ${seq_file} _ext.seq )_${cond}n${rep}.tsv
        all_meth_file_array=( ${all_meth_file_array[*]} ${meth_file} )
      done
    done

    cmd="
    methylation_metagene_script $( echo ${seq_files_all[*]} | tr ' ' ',' ) $( echo ${all_meth_file_array[*]} | tr ' ' ',' ) ${out_meth_meta_dir} ${SS} ${metagene_param} ${prefixes_list} ${expSeq_cond} ${expSeq_rep} ${comp_pair} ${min_reads} \
      --nbr-cores $NBRCORES --fig-prefix ${fig_prefix} ${color_pallette}
    "
    echo "$cmd"
    eval "$cmd"

  done

}
