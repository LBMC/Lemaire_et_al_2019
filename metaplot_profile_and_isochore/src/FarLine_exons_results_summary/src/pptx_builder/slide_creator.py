#!/usr/bin/python3

slide_layout = prs.slide_layouts[ 6 ]
slide = prs.slides.add_slide( slide_layout )

# textbox = slide.shapes.add_textbox( left=Inches( 0 ), top=Inches( 0 ), width=Inches( 6 ), height=Inches( 0.2 ) )
# textbox.text = 'exon list: ' + exon_list_name
# textbox.text_frame.paragraphs[0].runs[0].font.size = font_size_small

shapes = slide.shapes
shape = shapes.add_shape(
    MSO_SHAPE.RECTANGLE, left=slide_marge_left, top=Inches( 0.15 ), width=Inches( 6 ), height=Inches( 0.2 )
)
shape.fill.background()
shape.text = exon_list_name
shape.text_frame.paragraphs[0].alignment = PP_ALIGN.LEFT
shape.text_frame.paragraphs[0].font.size = Pt( 10 )
shape.text_frame.paragraphs[0].font.color.rgb = RGBColor( 0, 0, 0 )

textbox = slide.shapes.add_textbox( left=current_x, top=current_y, width=Inches( 1.8 ), height=height_tbx )
p = textbox.text_frame.paragraphs[0]
p.text = page_title
p.alignment = PP_ALIGN.LEFT
p.runs[0].font.size = font_size_big
