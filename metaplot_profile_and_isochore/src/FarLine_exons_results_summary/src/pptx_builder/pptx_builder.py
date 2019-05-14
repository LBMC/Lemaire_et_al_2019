#!/usr/bin/python3

from pptx import Presentation
from pptx.util import Inches,Pt
from pptx.enum.text import PP_ALIGN
from pptx.enum.shapes import MSO_SHAPE
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_THEME_COLOR
import sys
from os import makedirs,path,walk
from subprocess import getoutput

## define a text of commands (pseudo-functionalize execution of a file)
exec_text = '''with open( file_path, 'r' ) as fich:
    fich_text = ''.join( [ xxx for xxx in fich ] )
    pass
exec( fich_text )
'''

## define directory variables
script_dir = path.dirname( path.realpath( __file__ ) )
fig_dir = sys.argv[ 1 ] #~ '/mnt/d/home/inhiria/work_temp/id_card_temp/results/figures'
pptx_name = sys.argv[ 2 ] #~ './pptx_test'
exon_list_name = sys.argv[ 3 ] #~ 'siDNMT3b-siGL2'
exon_list_dir = sys.argv[ 4 ]
conf_file = sys.argv[ 5 ]

## load the configuration file
file_path = conf_file
exec( exec_text )

# print( in_exon_list )
in_exon_list = [ str( xxx ) for xxx in in_exon_list ]
out_off_list = [ str( xxx ) for xxx in out_off_list ]
cov_in_exon = str( cov_in_exon )
cov_out_off = str( cov_out_off )

## initialize pptx
prs = Presentation()
prs.slide_width = Inches( 8.27 )
prs.slide_height = Inches( 11.69 )

slide_marge_left = Inches( 0.39 )
slide_marge_right = slide_marge_left
slide_marge_top = slide_marge_left + Inches( 0.1 )
slide_marge_bottom = slide_marge_left

height_tbx = Inches( 0.3 )
img_height = Inches( 2.45 )
font_size_big = Pt( 16 )
font_size_small = Pt( 12 )

## exon_feature page
# try:
#     file_path = script_dir + '/exon_feature_page.py'
#     exec( exec_text )
# except:
#     print( '--- no exon_feature page', file=sys.stderr )
file_path = script_dir + '/exon_feature_page.py'
exec( exec_text )

## base_fraction page
try:
    file_path = script_dir + '/base_fraction_page.py'
    exec( exec_text )
except:
    print( '--- no base_fraction page', file=sys.stderr )

## siGL2_signal page
try:
    file_path = script_dir + '/siGL2_signal_page.py'
    exec( exec_text )
except:
    print( '--- no siGL2_signal page', file=sys.stderr )

## siPP-siGL2_signal page
try:
    file_path = script_dir + '/siPP-siGL2_signal_page.py'
    exec( exec_text )
except:
    print( '--- no siPP-siGL2 signal page', file=sys.stderr )

## public signal page
try:
    file_path = script_dir + '/public_signal_page.py'
    exec( exec_text )
except:
    print( '--- no public signal page', file=sys.stderr )

# ## delta_Psi page
# try:
#     file_path = script_dir + '/delta_Psi_page.py'
#     exec( exec_text )
# except:
#     print( '--- no delta_Psi page', file=sys.stderr )


out_file = pptx_name + '.pptx'
print( out_file )
prs.save( out_file )
