import svgwrite

import stem
from stem import Stem
from loop import Loop
from canvas import Canvas
from util import get_stems, flatten, consecutive_segments, color, seq_for

# General theory:
#
# We are interested in the relative arrangement of a set of stems.
#
# This is the most important problem to solve -- the placement of loop nucleotides
# is both less geometrically constrained by "traditional publication-quality aesthetics"
# and easier to resolve (in the sense of a loop-modeling problem, really).
#
# This is why we're first trying to *develop an API* against which to solve this problem.
#
# Stems have a position and orientation relative to the 'root stem' (at the moment, just
# stem 0; at the moment, the position and orientation *cannot* vary). So we can compute
# 'stem internal coordinates' -- trivial -- then transform them.

def apicals( stems, loops ):
	"""
	Finds all the apical loops in the given set of loops. All apical loops stretch
	from i+1 to j-1 for some base pair in a stem.
	"""
	
	# AMW: for draw_apical_loop_reasonably_with_respect_to_stem, we actually
	# want + need the stem_idx for now.
	apicals = []
	for loop in loops:
		print "Apical?", loop
		for stem_idx, stem in enumerate(stems):
			print "\t",stem.numbers
			if (min(loop)-1,max(loop)+1) in stem.numbers:
				apicals.append([(stem_idx, stem), loop])
				break

	return apicals

def junctions( stems, loops ):
	"""
	Junctions aren't just 'non-apical'; they're also 'non-tail'
	So instead of just taking the non apical loops, we instead have to ensure that our 
	candidate goes from max(stem1)+1 to min(stem2)-1
	"""
	junctions = []
	for loop in loops:
		for stem1_idx, stem1 in enumerate(stems):
			for stem2_idx, stem2 in enumerate(stems):
				if stem1_idx == stem2_idx: continue
				if stem2_idx != stem1_idx + 1: continue
				if max(flatten(stem1.numbers))+1 == min(loop) and min(flatten(stem2.numbers))-1 == max(loop):
					junctions.append([(stem1_idx, stem1), (stem2_idx, stem2), loop])

	print junctions
	return junctions

def draw_secstruct(stems, loops, seq, file_name='default.svg'):
	# init drawing
	dwg = svgwrite.Drawing(file_name)

	canvas = Canvas( dwg )
	for stem in stems: canvas.add_stem( stem )

	# Really, at that point we'll have to jointly figure out coordinates anyway. 
	for apical in [ loop for loop in loops if loop.apical ]: canvas.add_loop( apical )

	# Should the canvas be calling drawing functions
	# on the added stems/loops
	# call it 'render'
	for stem in stems: canvas.draw_stem( stem )

	for apical in [ loop for loop in loops if loop.apical ]:
		canvas.draw_apical_loop_reasonably_with_respect_to_stem( apical ) 

	# OK, by definition at the moment our y coords are zero 
	# and x coords are stem_idx * stem_offset_width + bp_offset_width and stem_idx * stem_offset_width 
	# so, interpolate
	for junction_loop in [ loop for loop in loops if not loop.apical ]:
		canvas.draw_junction_loop_reasonably_with_respect_to_stems( junction_loop )

	# save drawing
	dwg.save()
	return dwg

if __name__=="__main__":
	fn = 's_s_ss.svg'
	seq = 'ccccgcaaggggaucccauguucgcaug'
	ss  = '((((....))))....((((....))))'
    # Eventually: support chainbreaks
	stems = get_stems( ss, seq )
	# This produces old-style numerical stems. I want to turn them into
    # sequence stems. Use the sequence + numerical stem ctor.
	print stems
	stems = [ Stem( seq, stem ) for stem in stems ]

	print stems
	loops = []

	loop_nt = [ i+1 for i in xrange(len(seq)) if i+1 not in flatten([stem.numbers for stem in stems ])]
	#print "loop_nt", loop_nt
	loops = consecutive_segments( loop_nt )
	print loops
	real_loops = []
	for loop in loops:
		adjacent_stems = [ stem for stem in stems if min(loop)-1 in flatten(stem.numbers) or max(loop)+1 in flatten(stem.numbers) ]
		if len(adjacent_stems) == 0:
			print "loop", loop , "has no adjacent stems!"
			continue
		if len(adjacent_stems) == 1:
			real_loops.append( Loop( seq_for( loop, seq ), loop, adjacent_stems[0] ) )
		else: # junction. Don't cover 'tails' yet.
			real_loops.append( Loop( seq_for( loop, seq ), loop, adjacent_stems[0], False, adjacent_stems[1] ) )

	dwg = draw_secstruct(stems, real_loops, seq, file_name=fn)
