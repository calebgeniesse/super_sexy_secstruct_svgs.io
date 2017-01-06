import util
from util import flatten
from coordinate_frame import CoordinateFrame

# Thoughts:
# 1. A Stem or Loop shouldn't know details of how it is rendered.
# 2. An object containing multiple Stems and Loops might know how to orient them relative
# to each other
# 3. A Stem and Loop can carry information directly pertinent to rendering. For example, 
# a Stem doesn't just have a sequence (with a length); it has a number of base pairs
# with types.

WATSON_CRICK = 1

class Stem:
	def __init__( self, sequence, stem_pattern=None ):
		"""
		Initialize a stem with a given sequence. The length of the sequence 
		is used to size a variety of variables, such as its (relative/internal)
		coordinates.
		"""
		
		sequence = "".join( [c for c in sequence if c != ',' and c != ' ' ] )
		
		if stem_pattern is not None:
			sequence = "".join( [ sequence[i-1] for i in sorted(flatten(stem_pattern) ) ] )

		self.sequence = sequence 

		# AMW: Check later that this works for nontrivial examples
		# For example: will break on PK and intermolecular stems.
		#self.numbers  = sorted(flatten(stem_pattern) )
		# AMW: Actually, this may end up redoing the same organization as input
		# but I like being explicit
		numbers = sorted(flatten(stem_pattern) )
		self.numbers = [ ( numbers[i], numbers[len(numbers)-i-1] )  
		                    for i in xrange(len(numbers)/2 ) ]
		
		# Trim common spacers from sequence
		# Fold sequence into BPs
		self.base_pairs = [ ( sequence[i], sequence[len(sequence)-i-1] )  
		                    for i in xrange(len(sequence)/2 ) ]
		# All watson-crick for now.
		self.base_pair_types = [ WATSON_CRICK for x in self.base_pairs ]

		# First rotation notion. Pair with 1 (DOWN) or -1 (UP)
		self.coordinate_frame = CoordinateFrame( [ 0, 0 ], 1 )

