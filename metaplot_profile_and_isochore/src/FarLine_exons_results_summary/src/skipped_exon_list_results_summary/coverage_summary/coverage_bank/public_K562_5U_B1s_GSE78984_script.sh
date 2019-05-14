
  public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  expSeq_bw_dir=${data_dir}/bw_files/public_K562_5U_B1s_GSE78984
  source ${public_script_dir}/public_common_simpler.sh

  # maximum in mean metaplot
  mean_ymax=2.8
