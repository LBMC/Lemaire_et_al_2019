#!/usr/bin/python3

def unique_index ( liste ):
    #~ liste = [ 1,2,3,4,5,3,2,4,6,7 ]
    bouh = liste
    bouh.sort()
    bouh=[ bouh[ 0 ] ] + [ bouh[ xxx ] for xxx in list( range( 1, len( bouh ) ) ) if bouh[ xxx ] != bouh[ xxx - 1 ] ]
    return( [ liste.index( xxx ) for xxx in bouh ] )
