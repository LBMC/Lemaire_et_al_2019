#!/bin/bash

exon_coord_extraction_ASE_CE_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${exon_coord_extraction_ASE_CE_dir}/exon_coord_extraction.sh

feat_fasterdb=$1
ASE_list=$2
CE_list=$3

exon_coord_extraction "FasterDB" <( awk '( $11 == "alternative" )' $feat_fasterdb ) > $ASE_list
exon_coord_extraction "FasterDB" <( awk '( $11 == "constitutive" )' $feat_fasterdb ) > $CE_list
