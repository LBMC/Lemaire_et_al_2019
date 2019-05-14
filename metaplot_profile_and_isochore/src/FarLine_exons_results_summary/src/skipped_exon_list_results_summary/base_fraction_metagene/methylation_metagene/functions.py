#!/usr/bin/python3
import subprocess as sp
import sys

def load_file_whead( file_path, id_col=[], sep='\t' ):
    FL_tab = {}

    with open( file_path, 'r' ) as fich:
        junc_headers = fich.readline().strip().split( sep )
        #~ print( junc_headers )

        n = 0
        for line in fich:
            n += 1
            parts = line.strip().split( sep )
            junc_info = {}
            for head_idx, head in enumerate( junc_headers ):
                junc_info[ head ] = parts[ head_idx ]
                pass

            if id_col:
                FL_tab[ '_'.join( [ junc_info[ xxx ] for xxx in id_col ] ) ] = junc_info
            else:
                FL_tab[ n ] = junc_info
                pass

            pass

        pass

    # print( list( FL_tab.keys() )[ 0:10 ] )
    # print( FL_tab[ list( FL_tab.keys() )[ 0 ] ] )
    return( [ junc_headers, FL_tab ] )


def mpileup_fc( seq_dict, bam_file ):
    mpileup = sp.getoutput( ' '.join( [ '~/softwares/samtools-1.3.1/samtools', 'mpileup', '-q', '10', '-Q', '13', '-d', '100000', '-x', '-r', seq_dict[ 'coordinates' ], bam_file ] ) )
    mpileup_dict = {}
    for line in mpileup.split( '\n' ):
        parts = line.split( '\t' )
        try:
            mpileup_dict[ parts[ 0 ] + ':' + parts[ 1 ] ] = { 'chrom':parts[ 0 ], 'pos':int( parts[ 1 ] ), 'bases':parts[ 4 ] } #, 'total':int( parts[ 3 ] )
        except:
            continue
        pass

    return( mpileup_dict )


def CpG_sites_recover( seq_dict ):
    sequence = seq_dict[ 'sequence' ]
    chrom, borders = seq_dict[ 'coordinates'].split( ':' )
    strand = int( seq_dict[ 'strand' ] )

    # check strandness of the sequence to map to the reference sequence
    if strand == -1:
        ref_pos_idx = 1
    else:
        ref_pos_idx = 0
        pass
    ref_pos = int( borders.split( '-' )[ ref_pos_idx ] )

    # scan the sequence for the CpG sites
    CpG_list = []
    for pos in range( len( sequence ) - 1 ):
        if sequence[ pos:pos+2 ] == 'CG':
            CpG_list.append( ( chrom, ref_pos + int( ( strand - 1 ) / 2 ) + pos * strand, pos ) )
            pass
        pass

    return( CpG_list )


def CpG_meth_recover( CpG_list, mpileup_dict ):
    seq_meth_dict = {}
    for CpG_site in CpG_list:
        CpG_meth_dict = { 'C':{}, 'G':{} }
        # process the C of the CpG site
        C_coord = CpG_site[ 0 ] + ':' + str( CpG_site[ 1 ] )
        if C_coord in mpileup_dict:
            meth_dict = { 'raw':0, 'meth':0 }
            mpileup = mpileup_dict[ C_coord ]
            for xxx in range( len( mpileup[ 'bases' ] ) ):
                if mpileup[ 'bases' ][ xxx ] == 'C':
                    meth_dict[ 'meth' ] += 1
                elif mpileup[ 'bases' ][ xxx ] == 'T':
                    meth_dict[ 'raw' ] += 1
                    pass
                CpG_meth_dict[ 'C' ] = meth_dict
                pass
            pass

        # process the G of the CpG site
        G_coord = CpG_site[ 0 ] + ':' + str( CpG_site[ 1 ] + 1 )
        if G_coord in mpileup_dict:
            meth_dict = { 'raw':0, 'meth':0 }
            mpileup = mpileup_dict[ G_coord ]
            for xxx in range( len( mpileup[ 'bases' ] ) ):
                if mpileup[ 'bases' ][ xxx ] == 'g':
                    meth_dict[ 'meth' ] += 1
                elif mpileup[ 'bases' ][ xxx ] == 'a':
                    meth_dict[ 'raw' ] += 1
                    pass
                pass
            CpG_meth_dict[ 'G' ] = meth_dict
            pass

        seq_meth_dict[ C_coord ] = CpG_meth_dict
        pass

    return( seq_meth_dict )


def CpG_meth_out_write( seq_meth_dict ):
    # print header line
    print( '\t'.join( [ 'chrom', 'pos', 'meth_C', 'raw_C' ] ), file=sys.stdout )

    # print data lines
    for idx in seq_meth_dict:
        chrom, pos = idx.split( ':' )
        pos = int( pos )

        counts = seq_meth_dict[ idx ][ 'C' ]
        # print( idx )
        if ( 'meth' in counts ) and ( counts[ 'meth' ] + counts[ 'raw' ] ) != 0:
            print( '\t'.join( [ chrom, str( pos ), str( counts[ 'meth' ] ), str( counts[ 'raw' ] ) ] ), file=sys.stdout )
            pass

        counts = seq_meth_dict[ idx ][ 'G' ]
        pos = pos + 1
        if ( 'meth' in counts ) and ( counts[ 'meth' ] + counts[ 'raw' ] ) != 0:
            print( '\t'.join( [ chrom, str( pos ), str( counts[ 'meth' ] ), str( counts[ 'raw' ] ) ] ), file=sys.stdout )
            pass

        pass

    return( 0 )
