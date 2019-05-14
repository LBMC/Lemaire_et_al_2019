
## Script to source to build and set the general environment
root_dir=$1
src_dir=$2

## build the project struture
data_dir=${src_dir}/../data; ls -d ${data_dir} #~ mkdir -vp ${data_dir}
#~ genome_dir=${data_dir}/hg19/EnsEMBL_all; ls -d ${genome_dir} #~ mkdir -vp ${genome_dir}
results_ASE_CE_dir=${data_dir}/results_ASE-CE; ls -d ${results_ASE_CE_dir} > /dev/null
query_all_data_tab=${data_dir}/results_all_data_tab; ls -d ${query_all_data_tab} > /dev/null

res_dir=${root_dir}/results; mkdir -vp ${res_dir}

#~ ln -s ~/analyses_SEBASTIEN/data/all_faster_db_web.tsv ${data_dir}/all_faster_db_web.tsv 2> /dev/null
#~_ln -s ~/analyses_SEBASTIEN/data/tab_files/annotations/exons_genomiques_bis.tsv ${data_dir}/exons_genomiques_bis.tsv 2> /dev/null
#~_ln -s ~/analyses_SEBASTIEN/data/tab_files/annotations/genes.tsv ${data_dir}/genes.tsv 2> /dev/null
#~ for xxx in $( ls ~/analyses_SEBASTIEN/data/genome/hg19/EnsEMBL_all/*.fa ); do
#~   ln -s ${xxx} ${genome_dir}/$( basename ${xxx} ) 2> /dev/null
#~ done

## Set the environment
FEATURE_TAB=${data_dir}/all_faster_db_web.tsv
CODING_EXON_TAB=${data_dir}/hsapiens_exonsstatus_improved.csv
GENE_FDB_TAB=${data_dir}/genes.tsv
EXON_FDB_TAB=${data_dir}/exons_genomiques_bis.tsv
FASTA_DIR=${data_dir}/hg19/EnsEMBL_all
