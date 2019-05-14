#!/bin/bash

exon2SSregion () {
  local list_exon=$1
  local splicing_site=$2
  local in_exon=$3
  local out_off=$4

  echo -e "exon_id\tcoord\tstrand"

  if [[ "x$splicing_site" == x"3SS" ]]; then
    tail -n +2 ${list_exon} | \
      sed s/':'/'\t'/ | \
      sed s/'-'/'\t'/ | \
      awk -v in_exon=${in_exon} -v out_off=${out_off} '{
        OFS="\t";
        if ( $5 == "1" ) {
          print $1,$2":"$3-out_off"-"$3+in_exon-1,$5
        } else if ( $5 == "-1" ) {
          print $1,$2":"$4-in_exon+1"-"$4+out_off,$5
        }
      }'
  elif [[ "x$splicing_site" == x"5SS" ]]; then
    tail -n +2 ${list_exon} | \
      sed s/':'/'\t'/ | \
      sed s/'-'/'\t'/ | \
      awk -v in_exon=${in_exon} -v out_off=${out_off} '{
        OFS="\t";
        if ( $5 == "1" ) {
          print $1,$2":"$4-in_exon+1"-"$4+out_off,$5
        } else if ( $5 == "-1" ) {
          print $1,$2":"$3-out_off"-"$3+in_exon-1,$5
        }
      }'
  fi
}
