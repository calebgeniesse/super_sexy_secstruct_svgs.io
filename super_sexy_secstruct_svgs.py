import svgwrite

from rna_server_conversions import get_all_stems, join_sequence

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

stem_offset_width = 80
bp_offset_width = 30
bp_offset_height = 20

def color(nt):
	if nt == 'a': return 'yellow'
	if nt == 'c': return 'green'
	if nt == 'g': return 'red'
	if nt == 'u': return 'blue'
	
	print "WARNING: unrecognized nucleotide."

def draw_stems(stems, seq, dwg):
	# draw basepairs
	for stem_idx, stem in enumerate(stems):
		print  stem_idx, stem
		x_offset = stem_idx * stem_offset_width
		for bp_idx, [bp1, bp2] in enumerate(stem): 
			# need to compute some offet? ...
			y_offset = bp_idx * bp_offset_height
			dwg.add(dwg.line((x_offset + 10, y_offset), (x_offset + 20, y_offset), 
					stroke=svgwrite.rgb(10, 10, 16, '%')))
			dwg.add(dwg.text(seq[bp1-1], insert=(x_offset+0, y_offset), fill=color(seq[bp1-1])))#'blue'))
			dwg.add(dwg.text(seq[bp2-1], insert=(x_offset+bp_offset_width, y_offset), fill=color(seq[bp2-1])))#'yellow'))

			# numbers if multiple of 5
			# TODO: adjust size
			if bp1 % 5 == 0: dwg.add(dwg.text(bp1, insert=(x_offset+(-10), y_offset+5), fill=color(seq[bp1-1])))#'blue'))
			if bp2 % 5 == 0: dwg.add(dwg.text(bp2, insert=(x_offset+30+(10), y_offset+5), fill=color(seq[bp2-1])))#'blue'))


	return dwg

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
			print "\t",stem
			if [min(loop)-1,max(loop)+1] in stem:
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
				if max(flatten(stem1))+1 == min(loop) and min(flatten(stem2))-1 == max(loop):
					junctions.append([(stem1_idx, stem1), (stem2_idx, stem2), loop])

	print junctions
	return junctions

def draw_apical_loop_reasonably_with_respect_to_stem( apical_loop, stem_idx, stem, seq, dwg ): 
	"""
	In a better world, this stem would have attached coordinates.
	Instead I'm just passing its stem_idx and computing the offset.
	"""
	for loop_idx, loop_nt in enumerate(apical_loop):
		x_offset = stem_idx * stem_offset_width + loop_idx * 8
		# This would come out of stem length
		y_offset = len(stem) * bp_offset_height + 10
		dwg.add(dwg.text(seq[loop_nt-1], insert=(x_offset, y_offset), fill=color(seq[loop_nt-1])))#'blue'))
	return dwg

def draw_junction_loop_reasonably_with_respect_to_stems( junction_loop, stem1_idx, stem1, stem2_idx, stem2, seq, dwg ):
	#y1 = 0 - 10
	#y2 = 0 - 10
	x1 = stem1_idx * stem_offset_width + bp_offset_width
	x2 = stem2_idx * stem_offset_width
	for loop_idx, loop_nt in enumerate(junction_loop):
		print stem1_idx, stem2_idx, stem_offset_width,  bp_offset_width
		fraction_done_with_loop = float(loop_idx) / float(len( junction_loop ))
		print fraction_done_with_loop
		print (x1+(float(x2-x1)* fraction_done_with_loop ) )
		dwg.add(dwg.text(seq[loop_nt-1], insert=((x1+(float(x2-x1)* fraction_done_with_loop ) ), -10 ), fill=color(seq[loop_nt-1])))#'blue'))
	return dwg

def draw_secstruct(stems, loops, seq, file_name='default.out'):
	# init drawing
	dwg = svgwrite.Drawing(file_name)

	# draw struct here
	dwg = draw_stems(stems, seq, dwg)

	# Later, figure out whether particular loops are apical or junctions or whatever.
	# Really, at that point we'll have to jointly figure out coordinates anyway. 
	for (stem_idx, stem), apical_loop in apicals( stems, loops ):
		dwg = draw_apical_loop_reasonably_with_respect_to_stem( apical_loop, stem_idx, stem, seq, dwg ) 

	# We don't handle j1/2 yet
	# OK, by definition at the moment our y coords are zero 
	# and x coords are stem_idx * stem_offset_width + bp_offset_width and stem_idx * stem_offset_width 
	# so, interpolate
	for (stem1_idx, stem1), (stem2_idx, stem2), junction_loop in junctions( stems, loops ):
		dwg = draw_junction_loop_reasonably_with_respect_to_stems( junction_loop, stem1_idx, stem1, stem2_idx, stem2, seq, dwg ) 

	# save drawing
	dwg.save()
	return dwg

if __name__=="__main__":
	fn = 's_s_ss.svg'
	seq = 'ccccgcaaggggaucccauguucgcaug'
	ss  = '((((....))))....((((....))))'
	stems = get_all_stems( ss, [], seq )
	print stems
	loops = []
	# loops are consecutive sequences of non-stem residues
	def flatten( obj ):
		arr = []
		try:
			for i in obj:
				arr.extend( flatten( i ) )
			return arr
		except:
			return obj

	def consecutive_segments( lst ):
		lists = []
		newlist = []
		for idx, item in enumerate(lst):
			#print lst, lists, newlist
			if idx == 0:
				newlist.append(item)
				continue
			elif item - lst[idx-1] > 1:
				lists.append(newlist)
				newlist = []
				newlist.append(item)
				continue
			else: newlist.append(item)
		lists.append(newlist)
		
		if lists[-1] == []: return lists[:-1]
		return lists

	loop_nt = [ i+1 for i in xrange(len(seq)) if i+1 not in flatten(stems) ]
	#print "loop_nt", loop_nt
	loops = consecutive_segments( loop_nt )
	print loops
	dwg = draw_secstruct(stems, loops, seq, file_name=fn)
