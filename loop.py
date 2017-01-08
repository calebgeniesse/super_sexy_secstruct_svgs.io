import util
from util import flatten
from coordinate_frame import CoordinateFrame
from nucleotide import Nucleotide

class Loop:
	def __init__( self, sequence, numbers, stem1, apical=True, tail=False, stem2=None ):
		"""
		Initialize a loop with a particular sequence. Loops are associated with stems
		(we will not create SS diagrams for bare loops); either one or two.
		"""
		self.sequence = sequence
		self.numbers = numbers

		self.nucleotides = { numbers[i]: Nucleotide( sequence[i], numbers[i] ) for i in xrange(len(sequence)) }

		self.stem1 = stem1
		self.apical = apical
		self.tail = tail
		self.stem2 = stem2
		# The coordinate frame of a loop is intrinsically more complicated, since
		# a loop is tracing out a path and even if stems can only be dull things,
		# we still have notions of beginning and end here.
		#
		# This will also be good practice for manipulating more complicated 
		# coordinate frames later.
		# Note: orientation ignored.
		self.coordinate_frame = CoordinateFrame( [0,0], 1, [0,0] )
	
