#!/bin/bash

exons_list2BED6 () {
  tail -n +2 ${1:-/dev/stdin} | sed s/":"/'\t'/ | sed s/"-"/'\t'/ | awk '{ OFS="\t"; print $2,$3-1,$4,$1,0,$5 }' | sort -k1,1 -k2,2n
}
