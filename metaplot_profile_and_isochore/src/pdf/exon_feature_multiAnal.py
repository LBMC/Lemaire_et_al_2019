import glob

def img_pane_fc( image_file, slide, left_x, top_x, height ):
    try:
        img = slide.shapes.add_picture( image_file, left=left_x, top=top_x, height=height )
        # legend_hide( img, slide )
        obj_width = img.width
    except:
        textbox = slide.shapes.add_textbox( left=left_x, top=top_x, width=Inches( 0.3 ), height=height_tbx )
        textbox.text = 'no image'
        print( '---Error: missing `' + image_file + '`', file=sys.stderr )
        obj_width = siGL2_signal_img_height
        pass
    pass

def col_names_fc( slide, slide_marge_left, siGL2_signal_img_height, current_y, res_dir_list, offset_x ):
    current_y = current_y + height_tbx
    current_x = slide_marge_left + siGL2_signal_img_height / 2
    for res_dir_idx, res_dir in enumerate( res_dir_list ):
        current_pane_x = current_x + res_dir_idx * ( offset_x + siGL2_signal_img_height * 2 )
        textbox = slide.shapes.add_textbox( left=current_pane_x, top=current_y, width=siGL2_signal_img_height, height=height_tbx )
        textbox.text = exon_pos_list[ res_dir_idx ]
        textbox.text_frame.paragraphs[0].runs[0].font.size = font_size_small
        pass

    return( ( current_x, current_y ) )


# initialize position of the cursor
current_y = slide_marge_top
current_x = slide_marge_left

# create the slide
slide_layout = prs.slide_layouts[ 6 ]
slide = prs.slides.add_slide( slide_layout )
new_slide = True

# add title
textbox = slide.shapes.add_textbox( left=current_x + Inches( 1 ), top=current_y, width=prs.slide_width - Inches( 2 ), height=height_tbx )
p = textbox.text_frame.paragraphs[0]
p.text = page_title
p.alignment = PP_ALIGN.CENTER
p.font.size = font_size_big

# add the figures
inter_img_space_y = 0 #Inches( 0.1 )
offset_x = Inches( 0.3 )
siGL2_signal_img_height = floor( 1 * img_height / 2 )
second_column_x = prs.slide_width - slide_marge_right - 2 * siGL2_signal_img_height

fig_part_y = current_y + height_tbx + Inches( 0.1 ) # Inches( 1 )


#   # first column
#   #   # title of the column
current_y = fig_part_y #+ Inches( 0.1 )

#   #   # exon number
current_x, current_y = col_names_fc( slide, slide_marge_left, siGL2_signal_img_height, current_y, res_dir_list, offset_x )

#   #   # image panes
nb_meas = -1
for measure_idx, measure in enumerate( measure_list ):
    nb_meas += 1
    if nb_meas % 4 != 0 and not new_slide:
        current_y += siGL2_signal_img_height + 1/2 * height_tbx
    else:
        current_y = current_y + height_tbx #2 * height_tbx
        new_slide = False
        pass
    current_x = slide_marge_left

    textbox = slide.shapes.add_textbox( left=current_x + siGL2_signal_img_height / 2, top=current_y, width=siGL2_signal_img_height, height=height_tbx )
    textbox.text = measure
    textbox.text_frame.paragraphs[0].runs[0].font.size = font_size_small

    current_y += height_tbx
    for res_dir_idx, res_dir in enumerate( res_dir_list ):
        current_pane_x = current_x + res_dir_idx * ( offset_x + siGL2_signal_img_height * 2 )
        glob_pattern = res_dir + '/' + feature_dir + '/*_' + measure + '_boxplot.png'
        print( glob_pattern )
        image_file = glob.glob( glob_pattern )[ 0 ]
        img_pane_fc( image_file, slide, current_pane_x, current_y, siGL2_signal_img_height )
        pass

    if nb_meas == 3:
        # initialize position of the cursor
        current_y = slide_marge_top
        current_x = slide_marge_left

        # create the slide
        slide_layout = prs.slide_layouts[ 6 ]
        slide = prs.slides.add_slide( slide_layout )
        new_slide = True

        # add title
        textbox = slide.shapes.add_textbox( left=current_x + Inches( 1 ), top=current_y, width=prs.slide_width - Inches( 2 ), height=height_tbx )
        p = textbox.text_frame.paragraphs[0]
        p.text = page_title
        p.alignment = PP_ALIGN.CENTER
        p.font.size = font_size_big

        # add the figures
        inter_img_space_y = 0 #Inches( 0.1 )
        offset_x = Inches( 0.3 )
        siGL2_signal_img_height = floor( 1 * img_height / 2 )
        second_column_x = prs.slide_width - slide_marge_right - 2 * siGL2_signal_img_height

        fig_part_y = current_y + height_tbx + Inches( 0.1 ) # Inches( 1 )

        current_y = fig_part_y #+ Inches( 0.1 )
        current_x, current_y = col_names_fc( slide, slide_marge_left, siGL2_signal_img_height, current_y, res_dir_list, offset_x )
        pass
    pass
