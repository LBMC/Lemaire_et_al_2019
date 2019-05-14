#!/bin/bash
public_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  
  expSeq_bw_dir=${data_dir}/bw_files/RNA-Seq_Ctrl_SplicingLore
  encode_sign=HeLA_U2AF65_ctrl_RNASeq_siNC
  source ${public_script_dir}/public_common_encode.sh
  
