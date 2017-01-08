import svgwrite

import stem
from stem import Stem
from loop import Loop
from canvas import Canvas
from util import get_stems, flatten, consecutive_segments, color, seq_for
import numpy as np
import math
import copy

# test
#from util import loop_interpolate
#x1, y1 = [-1,0]
#x2, y2 = [ 1,0]
#print loop_interpolate( x1,y1,x2,y2, 0.75, 1.0/9.0 )
#exit()

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

def get_loops( seq, stems ):
	loop_nt = [ i+1 for i in xrange(len(seq)) if i+1 not in flatten([stem.numbers for stem in stems ])]
	#print "loop_nt", loop_nt
	loops = consecutive_segments( loop_nt )
	print loops
	real_loops = []
	for loop in loops:
		adjacent_stems = [ stem for stem in stems if min(loop)-1 in flatten(stem.numbers) or max(loop)+1 in flatten(stem.numbers) ]
		apical_stems = [ stem for stem in stems if min(loop)-1 in flatten(stem.numbers) and max(loop)+1 in flatten(stem.numbers) ]
		if len(adjacent_stems) == 0:
			print "loop", loop , "has no adjacent stems!"
			continue
		if len(apical_stems) == 1:
			real_loops.append( Loop( seq_for( loop, seq ), loop, apical_stems[0] ) )
		elif len(adjacent_stems) == 1: # tail
			real_loops.append( Loop( seq_for( loop, seq ), loop, adjacent_stems[0], False, True ) )
		else: # junction. Don't cover 'tails' yet.
			real_loops.append( Loop( seq_for( loop, seq ), loop, adjacent_stems[0], False, False, adjacent_stems[1] ) )
	return real_loops


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


# constants used below. May end up functions of font.
NT_DISTANCE = 20
# Imagine having different spring constants for stem vs loop vs PK
spring_constant = 1

def score( canvas ):
	"""
	Assign a score to a configuration of nucleotides.
	"""

	def distance( nt1, nt2 ):
		return math.sqrt( ( nt1.x - nt2.x ) ** 2 + ( nt1.y - nt2.y ) ** 2 )

	def harmonic_penalty( dist, ideal, spring_constant ):
		return spring_constant * ( dist - ideal ) ** 2

	score = 0

	# 1. Nucleotides should be a particular distance from what's adjacent to them in sequence.
	for seqpos in canvas.nucleotides.keys():
		if seqpos + 1 in canvas.nucleotides.keys():
			score += harmonic_penalty( distance( canvas.nucleotides[seqpos], canvas.nucleotides[seqpos+1] ), NT_DISTANCE, spring_constant )

	# 2. Nucleotides should be at least 15 from ALL other nts. Don't double-count.
	for seqpos1 in canvas.nucleotides.keys():
		for seqpos2 in canvas.nucleotides.keys():
			if seqpos1 >= seqpos2: continue
			if distance( canvas.nucleotides[seqpos1], canvas.nucleotides[seqpos2] ) < 15: score += 100
			

	return score

def perturb( canvas ):
	"""
	Random perturbation to each NT. Gaussian?
	Probably copying issues because Python.
	"""

	new_canvas = copy.deepcopy(canvas)#Canvas(canvas)
	print new_canvas.nucleotides.keys()
	for seqpos in new_canvas.nucleotides.keys():
		#print "Moving nt %d from ( %f, %f )... " % ( seqpos, new_canvas.nucleotides[seqpos].x, new_canvas.nucleotides[seqpos].y ),
		new_canvas.nucleotides[seqpos].x += np.random.normal(0,5)
		new_canvas.nucleotides[seqpos].y += np.random.normal(0,5)
		#print "... to ( %f, %f )!" % ( new_canvas.nucleotides[seqpos].x, new_canvas.nucleotides[seqpos].y )

	return new_canvas

def metropolis( new, old, temp ): 
	return np.random.uniform() < math.exp( ( old - new ) / temp )

def mc( canvas ):
	"""
	This is an in-progress framework for doing MCMC simulations on a canvas.
	The idea is that you score NT configurations, update, etc.
	"""
	cycles = 100
	for x in xrange(cycles):
		old_score = score(canvas)
		new_canvas = perturb(canvas)
		new_score = score(new_canvas)
		if new_score < old_score or metropolis( new_score, old_score, temp=1 ):
			# accept
			print "Accepted perturbation from %f to %f." % ( old_score, new_score )
			canvas = new_canvas

		else:
			print "Rejected perturbation from %f to %f." % ( old_score, new_score )
	return canvas

if __name__=="__main__":
	fn = 's_s_ss.svg'
	seq = 'ccaaccgcaagguuggaucccauguuaaaaacgcaug'
	ss  = '((((((....))))))....((((.........))))'
	#seq = 'ccaaccgcaagguuggaucccauguucgcaug'
	#ss  = '((((((....))))))....((((....))))'
	#seq = 'ccccgcaaggggaucccauguucgcaug'
	#ss  = '((((....))))....((((....))))'
    # Eventually: support chainbreaks

	# This produces old-style numerical stems. I want to turn them into
    # sequence stems. Use the sequence + numerical stem ctor.
	stems = [ Stem( seq, stem ) for stem in get_stems( ss, seq ) ]

	# Find loops by process of elimination
	loops = get_loops( seq, stems )

	dwg = svgwrite.Drawing( fn )

	canvas = Canvas( dwg )
	for stem in stems: canvas.add_stem( stem )

	
	# Inform the canvas that stem 1 and stem 2 are coaxial.
	# This doesn't make that great semantic sense, but the canvas
	# is the only stem-collection type object so it is aware of this
	# sort of thing.

	# Also, how is the API user really going to know?
	# You'd imagine they'd actually want to set it based on
	# labels containing residues or something...
	#canvas.set_stems_coaxial( 0, 1 )



	# We have to do this before adding any loops because (initial) positions
	# are calculated after addition time.
	for loop in loops: canvas.add_loop( loop )
	
	canvas = mc(canvas)
	
	canvas.render()
