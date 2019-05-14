#!/bin/bash

FasterDB_trans () {
  tail -n +2 ${1:-/dev/stdin} | awk '{ OFS="\t"; print $1"_"$2,$4 }' | sed s/'\tchr'/'\t'/
}

exon_coord_extraction () {
  local format=$1 #~ FarLine or FasterDB
  local liste_exons=$2
  local bed_exons_genomiques_bis=$3 #~ should be a file named exons_genomiques_bis.bed


  local sep="$( printf '\t' )"
  local champ1=1
  local champ2=1

  if [[ "x$format" == x"FarLine" ]]; then
    echo -e "exon_id\tcoord\tstrand"
    join -t "$sep" --header -1 $champ1 -2 $champ2 \
      <(
        tail -n +2 $liste_exons | awk '{ OFS="\t"; print $2,$4 }' | sort -t "$sep" -k ${champ1},${champ1}
       ) \
      <(
        tail -n +2 $bed_exons_genomiques_bis | awk '{ OFS="\t"; print $1,$10 }' | sort -t "$sep" -k ${champ2},${champ2}
       )

  elif [[ "x$format" == x"FasterDB" ]]; then
    echo -e "exon_gene_pos\tcoord\tstrand"

    local all_text="$( more $liste_exons )"
    local id="$( echo "$all_text" | FasterDB_trans | cut -f 1 -d "$sep" )"
    local coord="$( echo "$all_text" | FasterDB_trans | cut -f 2 -d "$sep" )"
    local coord_sep=':'

    paste \
      <( echo "$id" ) \
      <( paste -d ':' \
        <( echo "$coord" | cut -f 1 -d "$coord_sep" ) \
        <( echo "$coord" | cut -f 2 -d "$coord_sep" )
        ) \
      <( echo "$coord" | cut -f 3 -d "$coord_sep" | awk '{ print $0"1" }' )

  fi
}
