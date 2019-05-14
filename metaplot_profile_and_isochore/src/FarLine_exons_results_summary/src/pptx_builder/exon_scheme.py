#!/usr/bin/python3

## exon skipping scheme
init_left = Inches( 2.5 )
init_top = slide_marge_top
init_width = Inches( 1 )
init_height = Inches( 0.7 )
textbox_width = Inches( 0.4 )

shapes = slide.shapes
left = init_left
top = init_top
width = init_width
height = init_height

shape = shapes.add_shape(
    MSO_SHAPE.RECTANGLE, left, top, width, height
)
shape.text = 'exon\nbefore'
for xxx in range( len( shape.text_frame.paragraphs[0].runs ) ):
    shape.text_frame.paragraphs[0].runs[xxx].font.size = Pt( 14 )
    pass

textbox = shapes.add_textbox( left - textbox_width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '3\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )

textbox = shapes.add_textbox( left + width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '5\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )

left = left + width
top = top + height / 2
height = 0
shape = shapes.add_shape(
    MSO_SHAPE.LINE_INVERSE, left, top, width, height
)
shape.text = 'intron\nbefore'
for xxx in range( len( shape.text_frame.paragraphs[0].runs ) ):
    shape.text_frame.paragraphs[0].runs[xxx].font.size = Pt( 14 )
    shape.text_frame.paragraphs[0].runs[xxx].font.color.rgb = RGBColor( 0, 0, 0 )
    pass


left += width
top = init_top
height = init_height
shape = shapes.add_shape(
    MSO_SHAPE.RECTANGLE, left, top, width, height
)
shape.text = 'exon'
for xxx in range( len( shape.text_frame.paragraphs[0].runs ) ):
    shape.text_frame.paragraphs[0].runs[xxx].font.size = Pt( 14 )
    pass

textbox = shapes.add_textbox( left - textbox_width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '3\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )

textbox = shapes.add_textbox( left + width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '5\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )

left = left + width
top = top + height / 2
height = 0
shape = shapes.add_shape(
    MSO_SHAPE.LINE_INVERSE, left, top, width, height
)
shape.text = 'intron\nafter'
for xxx in range( len( shape.text_frame.paragraphs[0].runs ) ):
    shape.text_frame.paragraphs[0].runs[xxx].font.size = Pt( 14 )
    shape.text_frame.paragraphs[0].runs[xxx].font.color.rgb = RGBColor( 0, 0, 0 )
    pass


left += width
top = init_top
height = init_height
shape = shapes.add_shape(
    MSO_SHAPE.RECTANGLE, left, top, width, height
)
shape.text = 'exon\nafter'
for xxx in range( len( shape.text_frame.paragraphs[0].runs ) ):
    shape.text_frame.paragraphs[0].runs[xxx].font.size = Pt( 14 )
    pass

textbox = shapes.add_textbox( left - textbox_width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '3\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )

textbox = shapes.add_textbox( left + width, top + height, width = textbox_width, height = Inches( 0.3 ) )
textbox.text = '5\''
textbox.text_frame.paragraphs[0].runs[0].font.size = Pt( 14 )
