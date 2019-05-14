#!/bin/bash

exon_coord_extraction_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${exon_coord_extraction_script_dir}/exon_coord_extraction.sh

liste_exons_up=$1
liste_exons_down=$2
bed_exons_genomiques_bis=$3
exon_up_list=$4
exon_down_list=$5

exon_coord_extraction "FarLine" $liste_exons_up $bed_exons_genomiques_bis > $exons_up_list
exon_coord_extraction "FarLine" $liste_exons_down $bed_exons_genomiques_bis > $exons_down_list
