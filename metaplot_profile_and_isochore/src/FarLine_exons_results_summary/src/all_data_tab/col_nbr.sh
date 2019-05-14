#!/bin/bash

col_nbr () {
  # fichier=${1:-/dev/stdin}
  col_name=$1
  head -1 ${2:-/dev/stdin} | tr '\t' '\n' | grep -En ^"${col_name}"$ | cut -d ':' -f 1
}
