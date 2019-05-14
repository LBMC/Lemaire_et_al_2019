#!/usr/bin/python3

import sys

## Script taking a bed file of the exons of interest and extracting the corresponding all_data_tab lines
subGroup_file = open( sys.argv[ 1 ], 'r' )
all_data_tab_file = open( sys.argv[ 2 ], 'r' )

if "--header" in sys.argv:
    header_line = subGroup_file.readline()


all_data_array = {}

print( all_data_tab_file.readline(), end='' )
for line in all_data_tab_file:
    coords = line.split( '\t' )[ 3 ]
    all_data_array[ coords ] = line
    pass

coords_set = set()
for line in subGroup_file:
    coords = line.strip().split( '\t' )[ 0:3 ]
    coords = coords[ 0 ] + ':' + str( int( coords[ 1 ] ) + 1 ) + '-' + coords[ 2 ]
    if coords not in coords_set:
         coords_set.add( coords )
         print( all_data_array[ coords ], end='' )
         pass
