import math

def seq_for( loop, seq ):
	return "".join([ seq[idx-1] for idx in sorted(loop) ])

def color(nt):
	if nt == 'a': return 'yellow'
	if nt == 'c': return 'green'
	if nt == 'g': return 'red'
	if nt == 'u': return 'blue'
	
	print "WARNING: unrecognized nucleotide."

def interpolate( v1, v2, frac ): return v1 + float(v2-v1)*frac

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

def get_stems( line, sequence_for_fasta ):
	"""
	Not concerned with handling pseudoknots for now so let's just grab stems from our
	SS and sequence and get out of here. We are assuming that all stems are, in this
	sense, 'coequal.'

	A reasonable approach eventually would be to mandate that () is reserved for non-PK
	stems.
	"""

	chainbreak_pos = []
	spacers = [ ',', ' ' ]
	# Eventually: iterate
	left_bracket_char = '('
	right_bracket_char = ')'
	complement = {'a':['u', 'X[OMU]', 'X[PSU]', 'X[5MU]'], 'u':['a','g', 'X[OMG]'], 'c':['g'], 'g':['c','X[5MC]', 'X[OMC]', 'u', 'X[OMU]', 'X[PSU]', 'X[5MU]']};
	
	count = 0
	left_brackets = []
	pair_map = {}
	all_pairs = []
	
	# THIS, not the sequence itself, should be equal in size to the secstruct etc.
	fasta_entities = []
	i = 0
	while i < len(sequence_for_fasta):
		if sequence_for_fasta[i] in spacers: continue
		if i == len(sequence_for_fasta) - 1 or sequence_for_fasta[i+1] != '[':
			# simple case
			fasta_entities.append(sequence_for_fasta[i])
		else:
			entity = sequence_for_fasta[i]
			i += 1
			while sequence_for_fasta[i] != ']':
				entity += sequence_for_fasta[i]
				i += 1
			entity += ']'
			fasta_entities.append(entity)

		i += 1
	#print fasta_entities
				

	for i in range( len(line) ):
		if line[i] in spacers: continue
		count += 1
		if line[i] == left_bracket_char:  left_brackets.append( count )
		if line[i] == right_bracket_char:
			if len( left_brackets ) == 0:
				raise ValidationError( "Number of right brackets does not match left brackets" )
			res1 = left_brackets[-1]
			res2 = count
			del( left_brackets[-1] )
			pair_map[ res1 ] = res2
			pair_map[ res2 ] = res1
			all_pairs.append( [res1,res2] )
			if fasta_entities[ res2-1 ] in complement.keys() and len( fasta_entities ) > 0 and not ( fasta_entities[res1-1] in complement[ fasta_entities[res2-1] ] ):
				raise ValidationError( "Not complementary at positions %s%d and %s%d!"  % (fasta_entities[res1-1],res1,fasta_entities[res2-1],res2) )
	#print pair_map

	if ( len (left_brackets) > 0 ):
		raise ValidationError( "Number of right brackets does not match left brackets" )
	numres = count

	# Parse out stems
	already_in_stem = {}
	for i in range( numres ): already_in_stem[ i ] = 0

	stems = []
	stem_count = 0
	for i in range( 1, numres + 1 ):
		if not pair_map.has_key( i ) or already_in_stem[ i ]: continue

		k = i 
		stem_count += 1
		stem_res = []

		stem_res.append( [k, pair_map[k]] )
		already_in_stem[ k ] = 1
		already_in_stem[ pair_map[k] ] = 1

		while( pair_map.has_key( k + 1 ) and  pair_map[ k+1 ] == pair_map[ k ] - 1  and not already_in_stem[k+1] and (not k in chainbreak_pos) and ( not pair_map[k+1] in chainbreak_pos )  ):
			k += 1
			stem_res.append( [k, pair_map[k]] )
			already_in_stem[ k ] = 1
			already_in_stem[ pair_map[k] ] = 1

		stems.append( stem_res )
	return stems

def loop_interpolate( x1,y1,x2,y2, fraction_of_circle, fraction_done_with_loop, add_angle=0, loop_sense=1 ):
	"""
	Suppose (x1,y1) and (x2,y2) are two points separated by fraction_of_circle of
	arc. Then, return points fraction_done_with_loop along that arc between those
	two points.
	"""
	#    C_________D
	#    /\       /)|
	#   /  \     /) |   find coordinates of D
	#  /    \   )   |   given its fraction of total angle
	# /______\_/____|      away from B or whatever.
	#A    M   B     E
	#        x1,y1   xD,y(AB*mAE/mAB)

	# Nonsensically, loops are going in the wrong direction.
	# Equivalent to careful debugging of below is just taking the complement.

	fraction_done_with_loop = 1 - fraction_done_with_loop

	mMCB = (360-fraction_of_circle*360)/2
	mAB = math.sqrt( (x2-x1)**2 + (y2-y1)**2 )
	mMC = 0.5 * mAB * math.tan(math.radians(mMCB))
	mBC = 0.5 * mAB / math.sin(math.radians(mMCB))
	r = mBC

	# Thus all points satisfy x^2 + y^2 = mBC^2
	mBCD = fraction_done_with_loop * fraction_of_circle*360

	# MC perp MB, so C can be found: mMC in perp dir from av A, B
	unit = None
	if y1 == y2:
		unit = [0,1]
	elif x1 == x2:
		unit = [1,0]
	else:
		print y2-y1
		print x2-x1
		unit = [1,float(-1)/float((y2-y1)/float(x2-x1))]
	xC,yC = [ (unit[0] * mMC) + (x1+x2)/2, (unit[1] * mMC) + (y1+y2)/2]
	print xC, yC


	# xC + r * cos a
	# yC + r * sin a
	# a is angle from horizontal right
	# AMW TODO: for slightly nudged stems, this won't be -90
	# Add add_angle, which is 180 for flipped stems
	a = loop_sense*(-90 + add_angle + mMCB + mBCD)
	print "%0.2f %0.2f" % (xC + r * math.cos(math.radians(a)), yC + r * math.sin(math.radians(a)))
	return [ xC + r * math.cos(math.radians(a)), yC + r * math.sin(math.radians(a))]

