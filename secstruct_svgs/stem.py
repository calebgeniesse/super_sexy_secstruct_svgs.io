import util
from util import flatten
from coordinate_frame import CoordinateFrame
from base_pair import BasePair

# Thoughts:
# 1. A Stem or Loop shouldn't know details of how it is rendered.
# 2. An object containing multiple Stems and Loops might know how to orient them relative
# to each other
# 3. A Stem and Loop can carry information directly pertinent to rendering. For example, 
# a Stem doesn't just have a sequence (with a length); it has a number of base pairs
# with types.

WATSON_CRICK = 1

class Nucleotide:
	def __init__( self, one_letter_code, seqpos ):
		self.name = one_letter_code
		self.seqpos = seqpos
		self.x = 0
		self.y = 0
		self.bp_partner = None
		
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
		self.base_pairs = [ BasePair( 
								Nucleotide(sequence[i], numbers[i]), 
								Nucleotide(sequence[len(sequence)-i-1], numbers[len(sequence)-i-1]) )  
		                    for i in xrange(len(sequence)/2 ) ]
							
		# Build up nucleotides.
		self.nucleotides = {}
		for bp in self.base_pairs:
			self.nucleotides[ bp.nt1.seqpos ] = bp.nt1
			self.nucleotides[ bp.nt2.seqpos ] = bp.nt2

		# Right now, set internal coordinates. First nt of every stem will seem 
		# like a ref right now, but the canvas can connect things up as needed.
		# second half of seq ref is first half.
		for seqpos, nt in self.nucleotides.iteritems():
			if seqpos == min( numbers ): continue # min(numbers) stays None, 0, 0
			
			if self.numbers.index(seqpos) < len(self.numbers)/2:
				nt.ref_nt = self.nucleotides[ seqpos-1 ]
			else: #if self.numbers.index(seqpos) >= len(self.numbers)/2
				nt.ref_nt = nt.bp_partner
			nt.dx = nt.x - nt.ref_nt.x
			nt.dy = nt.y - nt.ref_nt.y
							
		# All watson-crick for now.
		self.base_pair_types = [ WATSON_CRICK for x in self.base_pairs ]

		# First rotation notion. Pair with 1 (DOWN) or -1 (UP)
		self.coordinate_frame = CoordinateFrame( [ 0, 0 ], 1 )
		print "in ctor, orientation ", self.coordinate_frame.orientation

