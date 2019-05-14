
  public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  expSeq_bw_dir=${data_dir}/bw_files/public_HeLaS3_B2s_GSE42951
  source ${public_script_dir}/public_common_simpler.sh

  # maximum in mean metaplot
  mean_ymax=7.3
