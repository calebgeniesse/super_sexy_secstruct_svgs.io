import super_sexy_secstruct_svgs
import sys
import svgwrite
from stem import Stem
from loop import Loop
from canvas import Canvas

def test_canvas_add_stem():
	dwg = svgwrite.Drawing( "foo" )
	canvas = Canvas( dwg )

	stem = Stem( 'acguacgu', [[1,8],[2,7],[3,6],[4,5]] )

	canvas.add_stem( stem )
	assert( len( canvas.nucleotides ) == 8 )

	
def test_canvas_rejects_stem_with_repeated_numbers():
	dwg = svgwrite.Drawing( "foo" )
	canvas = Canvas( dwg )

	stem = Stem( 'acguacgu', [[1,8],[2,7],[3,6],[4,5]] )
	stem2 = Stem( 'acguacgu', [[1,8],[2,7],[3,6],[4,5]] )

	with pytest.raises(Exception) as e_info:
		canvas.add_stem( stem )
		canvas.add_stem( stem2 )

def test_canvas_rejects_loop_with_repeated_numbers():
	dwg = svgwrite.Drawing( "foo" )
	canvas = Canvas( dwg )

	loop = Loop( 'acgu', [1,2,3,4], None )
	loop2 = Loop( 'acgu', [1,2,3,4], None )

	with pytest.raises(Exception) as e_info:
		canvas.add_loop( loop )
		canvas.add_loop( loop2 )

def test_sane_loop_order():

	# this stem is an impossibility of numbering -- just so that we can add to a canvas.
	stem = Stem( 'acguacgu', [[1,8],[2,7],[3,6],[4,5]] )


	loop = Loop( 'acgu', [1,2,3,4], stem )
	for loop_idx, loop_nt in enumerate( loop.nucleotides.keys() ):
		print loop_idx, loop_nt, loop.nucleotides[ loop_nt ].name

	dwg = svgwrite.Drawing( "foo" )
	canvas = Canvas( dwg )

	canvas.add_stem( stem )

	canvas.add_loop( loop )

	for loop_idx, loop_nt in enumerate( canvas.loops[0].nucleotides.keys() ):
		print loop_idx, loop_nt, loop.nucleotides[ loop_nt ].name

	
if __name__ == '__main__':
	test_sane_loop_order()