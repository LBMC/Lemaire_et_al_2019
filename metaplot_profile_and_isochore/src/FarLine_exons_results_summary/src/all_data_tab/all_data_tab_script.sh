#!/bin/bash

FarLine_analysis_dir=$1
res_dir=$2
exons_list_name=$3
GENE_FDB_TAB=$4
EXON_FDB_TAB=$5

## load scripts and set variables
all_data_tab_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${all_data_tab_script_dir}/plot2D_deltaPsi/plot2D_deltaPsi.sh

out_dir=${res_dir}/all_data_tab; mkdir -vp ${out_dir}


## query comparison part
# create all data table for the query comparison ($1:FarLine_analysis_dir) from tmp files
echo ">>> query all_data_tab"
query_all_data_tab=${out_dir}/all_data_${exons_list_name}.tab
  # check the all_data_tab is not provided
if [[ "x${@}" == x*"--query-all-data-tab"* ]]; then
  query_all_data_tab_input=$( echo "$@" | awk -F '--query-all-data-tab ' '{ print $2 }' | cut -d ' ' -f 1 )
  ln -fs ${query_all_data_tab_input} ${query_all_data_tab}
else
  Rscript ${all_data_tab_script_dir}/all_data_computer.r ${FarLine_analysis_dir} > ${query_all_data_tab}
fi


# create the tables of skipped exons, completing FarLine data with gene id and exon id from FasterDB, with up and down exons separated
echo ">>> separate Up and Down"
source ${all_data_tab_script_dir}/all_data_tab.sh ${out_dir} ${exons_list_name} ${query_all_data_tab} ${GENE_FDB_TAB} ${EXON_FDB_TAB}


## set of FarLine analyses part
echo ">>> FarLine analyses all_data_tab"
# load the FarLine analysis paths array
source ${all_data_tab_script_dir}/FarLine_analysis_paths.sh

# create all data table for each FarLine analysis
unset FarLine_all_data_array; declare -A FarLine_all_data_array
FL_all_data_dir=$( echo "$@" | awk -F '--FarLine-all-data-dir ' '{ print $2 }' | cut -d ' ' -f 1 )
for FarLine_analysis_name in ${!FarLine_analysis_path_array[@]}; do
  FarLine_all_data_array[${FarLine_analysis_name}]=${out_dir}/all_data_${FarLine_analysis_name}.tab
  if [[ -n "$FL_all_data_dir" ]]; then
    ln -fs ${FL_all_data_dir}/$( basename ${FarLine_all_data_array[${FarLine_analysis_name}]} ) ${FarLine_all_data_array[${FarLine_analysis_name}]}
  else
    Rscript ${all_data_tab_script_dir}/all_data_computer.r ${FarLine_analysis_path_array[${FarLine_analysis_name}]} > ${FarLine_all_data_array[${FarLine_analysis_name}]} 2> /dev/null
  fi
done


#~ # build the 2D plot of deltaPSI with the query Farline analysis
#~ echo ">>> plots deltaPsi"
#~ fig_name_table=''
#~ for FL_all_data_tab in ${!FarLine_all_data_array[@]}; do
#~   deltaPSI_plot_dir="${res_dir}/figures/deltaPSI_plot"
#~   fig_name_list="$( plot2D_deltaPsi ${query_all_data_tab} ${FarLine_all_data_array[$FL_all_data_tab]} ${exons_list_name} ${FL_all_data_tab} ${deltaPSI_plot_dir} )"
#~   fig_name_table="${fig_name_table}${FL_all_data_tab}\t${fig_name_list}\n"
#~ done
#~ echo -ne "${fig_name_table}" | sort > ${deltaPSI_plot_dir}/fig_name_table.tsv
