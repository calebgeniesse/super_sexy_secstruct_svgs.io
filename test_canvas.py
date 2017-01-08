import super_sexy_secstruct_svgs
import sys
import svgwrite
from stem import Stem
from canvas import Canvas

def test_canvas_add_stem():
	dwg = svgwrite.Drawing( "foo" )
	canvas = Canvas( dwg )

	stem = Stem( 'acguacgu', [[1,8],[2,7],[3,6],[4,5]] )

	canvas.add_stem( stem )
	assert( len( canvas.nucleotides ) == 8 )